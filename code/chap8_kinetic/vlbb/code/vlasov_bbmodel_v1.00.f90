!***********************************************************************************
!*  A Vlasov-Poisson/Vlasov-Ampere Solver, also for BB (Berk-Breizman) model.      *
!*          With and Without Krook Collison, Diffusion and Drag                    *
!*                vlasov_bbmodel.f90, 2011-08-26 07:15                             *
!*        Hua-sheng XIE, IFTS-ZJU, huashengxie@{gmail.com, zju.edu.cn}             *
!***********************************************************************************

module variable
	implicit real*8 (A-H,O-Z)
	integer :: stdout=0, method_efield=0, method_vlasov=0
	! stdout=0: output to display; 1: to file "diplay.txt".
	! method_efield: select method for solve Efield, i.e. how to solve the poisson
	!                or Ampere equation. see 'efield' subroutine.
	!                0. Fij99; 1. dE/dx; 2. Ampere; 3. 'tridag'.
	! method_vlasov: select method for solve Vlasov equation.

	real*8 :: gammad=0.d0                                ! damping rate for BB model
	integer,parameter :: m=63, n=64, nt=12000, nplot=400, ndiag=10
	! n: number of grid points for Xspace
	! m: number of grid points for v>0; Total Vspace grid points number is 2m+1
	! nt: number of time steps
	! nplot: number of time steps between f consecutive outputs.
	! ndiag: number of steps between diagnostic consecutive outputs.

	real*8,parameter :: rk0=.5d0, vmax=12.0d0, dt=0.001d0, pi=3.14159265358979323846
	! real*8 :: pi=2.d0*dacos(0.d0)
	! rk0: Wave number, in a periodic system gives plasma length.
	! vmax: maximum value of the velocity in Vthermal units
	! dt: time step (in plasma frequency units)
	real*8,parameter :: dv=vmax/dfloat(m), xl=4.d0*pi, dx=xl/dfloat(n)
	! xl: plasma length, length of the X box in Debye-length units
	! dv: mesh size in V space
	! dx: mesh size in X space

	! parameters for initial distribution function
	real*8,parameter :: alpha=0.01d0, fbeam=0.0
	! alpha: perturbation parameter
	! fbeam=1.0: two-stream/beam distribution; else: Maxwellian distribution

	real*8,parameter :: vtc=0.1d0, vtb=1.0d0, vdc=-0.0d0, vdb=3.0d0, nb=0.1d0
	! parameters for generating the initial distribution function, two drift
	! Maxwellian distributions, f0(v)=f0c(v)+f0b(v):
	!      f0c(v)=(1-nb)/(vtc*sqrt(2*pi))*exp(-((v-vdc)/(2*vtc))^2)
	!      f0b(v)=nb/(vtb*sqrt(2*pi))*exp(-((v-vdb)/(2*vtb))^2)

	real*8 :: dtodv=dt/dv
	real*8 :: dtodx=dt/dx
	! dtodv: dt diveded by dv. parameter used in the shift routine;=dt/dv.
	! dtodx: dt diveded by dx. parameter used in the shift routine;=dt/dx.
	integer :: it
	real*8 :: t, et0
	! it: index of time
	! t: actual time
	! et0: the initial total energy, used to calculate total energy change rate

	real*8,dimension(-m:m) :: v
	real*8,dimension(-1:n) :: x, ef, rho
	real*8,dimension(-1:n,-m:m) :: f
	! v: velocity values on the grid
	! x: value of X on the grid. The indices varies from 0 to N-1. The values -1 and
	!    n are used to simplify calculus at boundary.
	! ef: Electric field. Obtained after solving Poisson/Ampere equation.
	! rho: Density. Sum of f(i,j) over j.
	! f: distribution function values at the grid points.

	real*8,dimension(-m:m) :: sh, a
	integer,dimension(-m:m) :: is
	integer,dimension(-1:n,-m:m) :: ind
	! This four arrays precalculated in subroutine 'perp' to speed-up the shift in X
	! routine.

	real*8,dimension(0:n-1) :: Qm2, Sm2, Zm2, Qp4, Sp4, Zp4
	! This six arrays precalculated in subroutine 'perp' and used in 'efield0' subr

end module variable

!*---------------------------------------------------------------------------------*
program main

	use variable
	implicit real*8 (A-H,O-Z)

	real*8 ctime0, ctime1

	! obtain starting CPU time
	call CPU_TIME(ctime0)
	if(stdout==1) open(stdout,file="display.txt",status='replace')

	open(unit=15,file='history.out',form='formatted',status='replace')

	t=0.d0                                                            ! initial time

	call initial                                     ! introduces initial conditions
	call outputpara
	call efield                                                    ! initial E-Field

	! Time evolution
	do it=1,nt

		if(mod(it-1,ndiag)==0) call history
		if(mod(it-1,nplot)==0) then
			call plotf
			call plotef
		endif

		t=t+dt

		call vlasov                           ! solve the vlasov system for one step

	enddo
	
	call history

	call CPU_TIME(ctime1)
	write(stdout,*) 'Program CPU TIME=', ctime1-ctime0, ' s'
	if(stdout==1) close(stdout)
	close(15)

end program main

!*---------------------------------------------------------------------------------*
subroutine vlasov

	use variable
	implicit real*8 (A-H,O-Z)

	select case(method_vlasov)

		case(1)
			call shiftx1
			call efield
			call shiftv
			call shiftx1

		case default                                                         ! Fij99
			call shiftx     ! shift f(x,v,tn)=f(x-v dt/2,v,tn*) for half a time step
			call efield                        ! solution of Poisson/Ampere equation
			call shiftv          ! shift f(x,v,tn*)=f(x,v-E dt,tn**) for a time step
			call shiftx    ! shift f(x,v,tn)=f(x-v dt/2,v,tn**) for half a time step

	end select

end subroutine vlasov

!*---------------------------------------------------------------------------------*
subroutine initial

	use variable
	implicit real*8 (A-H,O-Z)
	real*8 :: cte

	! initial grid
	! ************************************
	do i=-1,n
		x(i)=dfloat(i)*dx                        ! initialization of the Xspace grid
	enddo
	do j=-m,m
		v(j)=dfloat(j)*dv						 ! initialization of the Vspace grid
	enddo
	! ************************************

	! initial distribution function f(x,v)
	! *************************************
	cte=1.d0/dsqrt(2.d0*pi)
	do i=0,n-1
		pert=alpha*dcos(rk0*x(i))                          ! pert: perturbation term
		! Homogeneous with a perturbation in X
		do j=-m,m
			if(fbeam==1) then
				! Two drift Maxwellian distributions, f0(v)=f0c(v)+f0b(v):
				!      f0c(v)=(1-nb)/(vtc*sqrt(2*pi))*exp(-((v-vdc)/(2*vtc))^2)
				!      f0b(v)=nb/(vtb*sqrt(2*pi))*exp(-((v-vdb)/(2*vtb))^2)
				f(i,j)=(cte*(1-nb)/vtc*dexp(-0.5d0*((v(j)-vdc)/(vtc))**2)+&
				 cte*nb/vtb*dexp(-0.5d0*((v(j)-vdb)/(vtb))**2))*(1.d0+pert)
			else
				! f(i,j): Maxwellian in V.
				f(i,j)=cte*dexp(-v(j)*v(j)/2.0d0)*(1.d0+pert)
			endif
		enddo
	enddo
	
	fsm=0.d0
	do i=0,n-1                                      ! sum on f(x,v) to ensure be 1.0
		do j=-m,m
			fsm=fsm+f(i,j)
		enddo
	enddo
	fsm=fsm*dv/dfloat(n)                                             ! sum of f(x,v)
	write(stdout,*) 'Sum of f0(v)=', fsm

	! Without the line below the result may become very inaccurate, especially when
	! f0 is not normalized well!!!
	f(:,:)=f(:,:)/fsm        ! ensure sum of f(x,v) to 1.0

	! *************************************

	! boundary conditions
	! *************************************
	do i=-1,n
		f(i,-m)=0.d0          ! distribution function last value in V is forced to 0
		f(i,m)=0.d0
	enddo
	do j=-m,m
		f(n,j)=f(0,j)         ! In X space in the present case we impose periodicity

		! For an open system an empty space is needed to aviod matter loss.
		f(-1,j)=f(n-1,j)
	enddo

	call preps                ! precalculations of arrays used in the shiftx routine
	call prepe               ! precalculations of arrays used in the efield0 routine

end subroutine initial

!!! Note: 'preps' haven't been fully understood yet, and should be improved
!*---------------------------------------------------------------------------------*
subroutine preps                   ! From Subroutine 'PREP' of Fijalkow1999([1]&[2])

	! prep: This subroutine generates at the begining of the program arrays, used 
	!       for the X shift.

	use variable
	implicit real*8 (A-H,O-Z)

	do j=-m,-1               ! Calculattion of the shift coefficients for negative V
		sh(j)=-v(j)*dtodx                           ! value of the shift, normalized
		ii=int(sh(j))                            ! the number of cells to be shifted
		sh(j)=sh(j)-dfloat(ii)           ! the shift to be perform in the given cell
		a(j)=-.25d0*(1.d0-sh(j))                             ! the shift coefficient
		do i=0,n
			ind(i,j)=mod(i+ii,n)    ! index of the cell where the shift is performed
		enddo
		ind(-1,j)=mod(ii-1+n,n)
		! Since the system is periodic the value is calculated modulo n. In non-
		! periodic case ind must be smaller than n.
	enddo

	do j=0,m                 ! Calculattion of the shift coefficients for positive V
		sh(j)=v(j)*dtodx                            ! value of the shift, normalized
		ii=int(sh(j))                            ! the number of cells to be shifted
		sh(j)=sh(j)-dfloat(ii)        !at this stage sh is the shift fractional part
		a(j)=.25d0*(1.d0-sh(j))                              ! the shift coefficient
		do i=0,n
			ind(i,j)=mod(i-ii,n)    ! index of the cell where the shift is performed
		enddo
		ind(-1,j)=mod(-ii-1+n,n)
	enddo
	do j=-m,-1
		is(j)=1                                            ! Index, positive for V<0
	enddo
	do j=0,m
		is(j)=-1                                          ! Index, negative for V>=0
	enddo

end subroutine preps

!!! Note: 'prepe' haven't been fully understood yet, and should be improved
!*---------------------------------------------------------------------------------*
subroutine prepe                   ! From Subroutine 'PREP' of Fijalkow1999([1]&[2])

	! prepe: This subroutine generates at the begining of the program arrays, used 
	!       for the EField0.

	use variable
	implicit real*8 (A-H,O-Z)

	a(0)=0.d0
	Qp4(0)=-.25d0      ! constants array needed to solved the spline of the potential
	Sp4(0)=-.25d0      ! in equation: E(i+1)+4*E(i)+E(i-1) = 3*(phi(i+1)-phi(i-1))/dx
	Zp4(n-1)=1d0      !                                   
	do i=1,n-1
		Qp4(i)=-1d0/(Qp4(i-1)+4d0)
		Sp4(i)=Sp4(i-1)*Qp4(i)
	enddo
	do i=n-2,0,-1
		Zp4(i)=Qp4(i)*Zp4(i+1)+Sp4(i)
	enddo

	Qm2(0)=.5d0       ! constants array needed to solved the spline of the potential
	Sm2(0)=.5d0       ! in equation: phi(i+1)-2*phi(i)+phi(i-1) =
	Zm2(n-1)=1d0      !                        dx^2*(rho(i-1)+10*rho(i)+rho(i+1))/12
	do i=1,n-1
		Qm2(i)=-1d0/(Qm2(i-1)-2d0)
		Sm2(i)=Sm2(i-1)*Qm2(i)
	enddo
	do i=n-2,0,-1
		Zm2(i)=Qm2(i)*Zm2(i+1)+Sm2(i)
	enddo
	
end subroutine prepe

!!! Note: 'shiftx' haven't been fully understood yet, and should be improved
!*---------------------------------------------------------------------------------*
subroutine shiftx                     ! Subroutine 'SHIFTX' of Fijalkow1999([1]&[2])

	! Subroutine 'shiftx' shift the distribution function f(x,v,t) for a half
	! a time step; f(x-v dt,v,tn)=f(x,v,tn*).
	! For this subroutine precomputated shift coefficients are used

	use variable
	implicit real*8 (A-H,O-Z)

	real*8,dimension(-1:n) :: y, df
	! y: Temporary storage for the distribution function to be shifted
	! df: Phase space "Fluid" shifted from cell i to i+1

	do j=-m,m
		do i=-1,n
			y(i)=f(ind(i,j),j)              ! circular (for periodic problems) shift
		enddo
		do i=0,n-1              ! computation of the transported phase space "Fluid"
			df(i)=sh(j)*(y(i)+a(j)*(y(i+1)-y(i-1)))
		enddo
		df(n)=df(0)        ! introdution of periodic boundary values for moved fluid
		df(-1)=df(n-1)
		do i=0,n-1       ! new distribution function (old value-mass loss+mass gain)
			f(i,j)=y(i)+df(i+is(j))-df(i)
		enddo
	enddo

	do j=-m,m                ! periodicity is maintainded at the boundary in X space
		f(n,j)=f(0,j)
		f(-1,j)=f(n-1,j)
	enddo

end subroutine shiftx

!!! Note: 'shiftv' haven't been fully understood yet, and should be improved
!*---------------------------------------------------------------------------------*
subroutine shiftv                     ! Subroutine 'SHIFTV' of Fijalkow1999([1]&[2])

	! Subroutine 'shiftv' shift the distribution function f(x,v,t) for a time 
	! step; f(x,v,tn*)=f(x,v-E dt,tn**).
	use variable
	implicit real*8 (A-H,O-Z)

	real*8,dimension(-m:m) ::  y, df
	! y: Temporary storage for the distribution function to be shifted
	! df: Phase space "mass" shifted from cell j to j+1

	do i=0,n-1                        ! loop on all the grid point in the variable X
		if(ef(i) .gt. 0.d0) then        ! shift for electric field positive values  
		                                ! previously calculated by a call to efield.
			shh=ef(i)*dtodv                         ! value of the shift, normalized
			jj=int(shh)         ! shift integer value. jj is the number of the cells
							    ! to be shifted. As the system is non periodic jj
							    ! must be smaller in absolute value than m.
			shh=shh-dfloat(jj)   ! The fractional part of the shift to be perform in
								 ! the given cell. The same for all velocities.
			aa=0.25d0*(1.d0-shh)                             ! the shift coefficient
			do j=-m,jj-1-m                            ! EOSHIFT of f by a factor -jj
				y(j)=0.d0                              ! Empty cells from the boader
			enddo

			do j=-m+jj,m-1
				y(j)=f(i,j-jj)                                         ! shift cells
			enddo

			do j=1-m,m-1
				df(j)=shh*(y(j)+aa*(y(j+1)-y(j-1)))! the shifted phase space "Fluid"
			enddo

			df(-m)=0.d0                         ! To avoid matter to quit the system
			df(m)=0.d0                                ! boundary values are set to 0
			y(m)=0.d0
			do j=1-m,m-1                                 ! new distribution function
				f(i,j)=y(j)+df(j-1)-df(j)        ! (old value-fluid loss+fluid gain)
			enddo

		else                              ! shift for electric field negative values
			shh=-ef(i)*dtodv                               ! shh: value of the shift
			jj=int(shh)                      ! the number of the cells to be shifted
			shh=shh-dfloat(jj)   ! the fractional part of the shift to be perform in
			                     ! the given cell. The same for all velocities.

			aa=-0.25d0*(1.d0-shh)                            ! the shift coefficient
			do j=1-m,m-jj                             ! EOSHIFT of f by a factor +jj
				y(j)=f(i,j+jj)                                       ! shifted cells
			enddo
			
			do j=m-jj+1,m
				y(j)=0.d0                              ! Empty cells from the border
			enddo
			y(-m)=0.d0
			do j=1-m,m-1
				df(j)=shh*(y(j)+aa*(y(j+1)-y(j-1)))! the shifted phase space "Fluid"
			enddo

			df(-m)=0.d0                         ! To avoid matter to quit the system
			df(m)=0.d0                                ! Boundary values are set to 0
			do j=1-m,m-1                                 ! new distribution function
				f(i,j)=y(j)+df(j+1)-df(j)        ! (old value-fluid loss+fluid gain)
			enddo

		endif
	enddo
	do i=-1,n                                   ! boundary in V space f must be null
		f(i,-m)=0.d0
		f(i,m)=0.d0
	enddo
	do j=-m,m                                     ! boundary in X space are periodic
		f(n,j)=f(0,j)
		f(-1,j)=f(n-1,j)
	enddo

end subroutine shiftv

!  Need to check
!*---------------------------------------------------------------------------------*
subroutine shiftx1

	use variable
	implicit real*8 (A-H,O-Z)

	do j=-m,m
		call upwindscheme(f(:,j),v(j),n,0.5d0*dtodx)   ! a half time step
	enddo

end subroutine shiftx1

!  Need to check
!*---------------------------------------------------------------------------------*
subroutine upwindscheme(u,a,n,dtodx)      ! Upwind Scheme for the advection equation
	
	! Solve the advection equation, du/dt+a*du/dx=0, upwind Scheme. The CFL
	! condition should be satisfied: c=|a*dt/dx|<=1.
	! Ref: http://en.wikipedia.org/wiki/Upwind_scheme
	! Performance not good: 1st-order, diffusive; 2nd&3rd, dispersive and errors
	! (caused by boundary condition ??).

	implicit real*8 (A-H,O-Z)
	real*8,dimension(-1:n) :: u, a, uup,udown, aup, adown

	! calculation of a_up and a_down
	do i=-1,n
		aup(i)=max(a(i),0.d0)
		adown(i)=min(a(i),0.d0)
	enddo

	! calculation of u_up and u_down
	select case (1)                                                   ! select order
		case(3)	                                                         ! 3rd-order
			udown(0)=(2.d0*u(1)+3.d0*u(0)-6.d0*u(-1)+u(n-2))/6.d0
			uup(0)=(-u(2)+6.d0*u(1)-3.d0*u(0)-2.d0*u(n-1))/6.d0
			udown(n-1)=(2.d0*u(0)+3.d0*u(n-1)-6.d0*u(n-2)+u(n-3))/6.d0
			uup(n-1)=(-u(1)+6.d0*u(n)-3.d0*u(n-1)-2.d0*u(n-2))/6.d0
			do i=1,n-2
				udown(i)=(2.d0*u(i+1)+3.d0*u(i)-6.d0*u(i-1)+u(i-2))/6.d0
				uup(i)=(-u(i+2)+6.d0*u(i+1)-3.d0*u(i)-2.d0*u(i-1))/6.d0
			enddo
		case(2)                                                          ! 2nd-order
			udown(0)=(3.d0*u(0)-4.d0*u(-1)+u(n-2))/2.d0
			uup(0)=(-u(2)+4.d0*u(1)-3.d0*u(0))/2.d0
			udown(n-1)=(3.d0*u(n-1)-4.d0*u(n-2)+u(n-3))/2.d0
			uup(n-1)=(-u(1)+4.d0*u(0)-3.d0*u(n-1))/2.d0
			do i=1,n-2
				udown(i)=(3.d0*u(i)-4.d0*u(i-1)+u(i-2))/2.d0
				uup(i)=(-u(i+2)+4.d0*u(i+1)-3.d0*u(i))/2.d0
			enddo
		case default 	                                                 ! 1st-order
			do i=0,n-1
				udown(i)=u(i)-u(i-1)
				uup(i)=u(i+1)-u(i)
			enddo

	end select

	! update u
	do i=0,n-1
		u(i)=u(i)-(aup(i)*udown(i)+adown(i)*uup(i))*dtodx
	enddo

	! periodical boundary condition in X
	u(-1)=u(n-1)
	u(n)=u(0)
	
end subroutine upwindscheme

!  Need to check
!*---------------------------------------------------------------------------------*
subroutine fluxbalance(u,a,n,dtodx) ! Flux balance scheme for the advection equation
	
	! Solve the advection equation, du/dt+a*du/dx=0, Flux balance method. The CFL
	! condition should be satisfied: c=|a*dt/dx|<=1.
	! Ref: Fij99. 2nd-order

	implicit real*8 (A-H,O-Z)
	real*8,dimension(-1:n) :: u, a, uup,udown, aup, adown

	! calculation of a_up and a_down
	do i=-1,n
		aup(i)=max(a(i),0.d0)
		adown(i)=min(a(i),0.d0)
	enddo

	! calculation of u_up and u_down,  2nd-order

	udown(0)=(3.d0*u(0)-4.d0*u(-1)+u(n-2))/2.d0
	uup(0)=(-u(2)+4.d0*u(1)-3.d0*u(0))/2.d0
	udown(n-1)=(3.d0*u(n-1)-4.d0*u(n-2)+u(n-3))/2.d0
	uup(n-1)=(-u(1)+4.d0*u(0)-3.d0*u(n-1))/2.d0
	do i=1,n-2
		udown(i)=(3.d0*u(i)-4.d0*u(i-1)+u(i-2))/2.d0
		uup(i)=(-u(i+2)+4.d0*u(i+1)-3.d0*u(i))/2.d0
	enddo

	! update u
	do i=0,n-1
		u(i)=u(i)-(aup(i)*udown(i)+adown(i)*uup(i))*dtodx
	enddo

	! periodical boundary condition in X
	u(-1)=u(n-1)
	u(n)=u(0)
	
end subroutine fluxbalance

!*---------------------------------------------------------------------------------*
subroutine efield                                  ! Select a method to solve Efield

	use variable
	implicit real*8 (A-H,O-Z)

	select case(method_efield)
		case(0)
			call efield0                                     ! Fijalkow1999([1]&[2])
		case(1)
			call efield1                           ! Direct finite difference method
		case(2)
			call efield2                                          ! Via Ampere's Law
		case(3)
			call efield3                                          ! 'tridiag' method
		case default
			call efield0
	end select

end subroutine efield

!*---------------------------------------------------------------------------------*
subroutine efieldcorr                                            ! Efield correction
	! Ensure Efield periodicity on boundaries and the mean field value equal to zero

	use variable
	implicit real*8 (A-H,O-Z)

	! Ensure the mean field value equal to 0
	! *************************************
	efmean=0.d0                                                   ! mean field value
	do i=0,n-1
		efmean=efmean+ef(i)  ! integration const. so the mean field value equal to 0
	enddo
	efmean=efmean/dfloat(n)
	do i=-1,n
		ef(i)=ef(i)-efmean                                    ! final electric field
	enddo

	ef(n)=ef(0)                                          ! periodicity on boundaries
	ef(-1)=ef(n-1)

end subroutine efieldcorr

! Note: To keep 'efield0' go well, all real variables should use double precision
! (10^-16) in the whole program.
!*---------------------------------------------------------------------------------*
subroutine efield0               ! Use subroutine 'EFIELD3' of Fijalkow1999([1]&[2])
	
	! should call 'prepe' in the main program before using this routine

	use variable
	implicit real*8 (A-H,O-Z)
	real*8 :: dens, densm, phin
	real*8,dimension(0:n-1) :: u, phi, w, d

	densm=0.d0
	do i=0,n-1                         ! f(x,v) sum on V to obtain electrons density
		dens=0.d0
		do j=-m,m
			dens=dens+f(i,j)
		enddo
		rho(i)=dens*dv                             ! rho: density at the grid points
		densm=densm+dens
	enddo
	densm=densm*dv/dfloat(n)           ! mean density, or the contribution from ions
	do i=0,n-1
		rho(i)=rho(i)-densm    ! minus ions part, to ensure the mean density be zero
	enddo

	dx2d12=dx**2/12.0d0
	d(0)=(rho(n-1)+1.0d1*rho(0)+rho(1))*dx2d12                    ! Smoothed density
	d(n-1)=(rho(n-2)+1.0d1*rho(n-1)+rho(0))*dx2d12
	do i=1,n-2
		d(i)=(rho(i+1)+1.0d1*rho(i)+rho(i-1))*dx2d12
	enddo
	u(0)=-0.5d0*d(0)                   ! Intermediate values to potential calculation
	do i=1,n-1                          ! By spline method (see G.Knorr et al, 1979)
		u(i)=-(d(i)-u(i-1))*Qm2(i)
	enddo

	w(n-1)=0.0d0
	do i=n-2,0,-1
		w(i)=Qm2(i)*w(i+1)+u(i)
	enddo

	! Here very sensitivity, one should use double precision (10^-16) for all real
	! variables in this program to keep the accuracy and avoid meeting '0 over 0'!!!
	phin=(d(n-1)-w(0)-w(n-2))/(Zm2(0)+Zm2(n-2)-2d0)
	
	phi(n-1)=phin   ! Potential values from spline

	do i=0,n-2
		phi(i)=Zm2(i)*phin+w(i)
	enddo

	tddx=3d0/dx
	do i=1,n-2
		d(i)=(phi(i+1)-phi(i-1))*tddx               ! RHS term for field calculation
	enddo
	d(0)=(phi(1)-phi(n-1))*tddx                     ! Field calculation Intermediate
	d(n-1)=(phi(0)-phi(n-2))*tddx                   ! values, again by spline
	u(0)=-.5d0*d(0)
	do i=1,n-1
		u(i)=-(d(i)-u(i-1))*Qp4(i)
	enddo
	w(n-1)=0d0
	do i=n-2,0,-1
		w(i)=Qp4(i)*w(i+1)+u(i)
	enddo
	efn=(d(n-1)-w(0)-w(n-2))/(Zp4(0)+Zp4(n-2)+4d0)
	ef(n-1)=efn
	do i=0,n-2
		ef(i)=Zp4(i)*efn+w(i)                    ! Electric field values from spline
	enddo

	! Ensure Efield periodicity on boundaries and the mean field value equal to zero
	call efieldcorr

end subroutine efield0

!*---------------------------------------------------------------------------------*
subroutine efield1             ! Use direct finite difference method to solve Efield

	! This is the simplest method to get Efield, considered ES1D, then, dE/dx=rho, 
	! which can be solved directly.
	! The equation is solved for the index i varying from 0 to n-1; Values of the
	! field at the points -1 and N are given here by periodicity. In non periodic 
	! cases these values must be imposed.

	use variable
	implicit real*8 (A-H,O-Z)

	do i=-1,n
		dens=0.d0
		do j=-m,m
			dens=dens+f(i,j)                                  ! sum of f(x,v) over V
		enddo
		rho(i)=dens*dv                           ! rho: density at the X grid points
	enddo

	ef(0)=0.d0
	do i=1,n-1                                 ! solve dE/dx=rho_elelctrons-rho_ions
		ef(i)=ef(i-1)+0.5d0*(rho(i-1)+rho(i)-2.0d0)*dx
	enddo

	! Ensure Efield periodicity on boundaries and the mean field value equal to zero
	call efieldcorr

end subroutine efield1

!!! Note: 'efield2' haven't been checked yet
!*---------------------------------------------------------------------------------*
subroutine efield2                   ! Solve the equation of Ampere's Law for Efield

	use variable
	implicit real*8 (A-H,O-Z)

	real*8,dimension(-m:m) :: f0
	real*8,dimension(-1:n) :: def  ! r.h.s of the Ampere's Law

	if(t==0.d0) call efield1                                ! Get the initial Efield

	! Calculate the spatially averaged distribution function
	! *************************************
	hh=0.d0
	do j=-m,m
		h=0.d0
		do i=0,n-1
			h=h+f(i,j)
		enddo
		f0(j)=h/dfloat(n)
		hh=hh+h
	enddo
!	f00=hh
	! *************************************
	
	! Calculate the Efield use the former time step value, via Ampere's Law
	! *************************************
	do i=0,n-1
		hh=0.d0
		do j=-m,m
		!	hh=hh-v(j)*(f(i,j)-f0(j))
			hh=hh-v(j)*f(i,j)
		enddo
		! def(i)=hh*dv
		def(i)=hh*dv
		ef(i)=ef(i)+(def(i)-gammad*ef(i))*dt
	enddo

	! Ensure Efield periodicity on boundaries and the mean field value equal to zero
	call efieldcorr

end subroutine efield2

!*---------------------------------------------------------------------------------*
subroutine efield3   ! Use 'tridag' from Numerical Recipes to solve Phi, then Efield

	use variable
	implicit real*8 (A-H,O-Z)

	real*8,dimension(-1:n) :: d, phi, gam
	real*8 :: aa, bb, cc, bet

	! Calculation for density
	! *************************************
	densm=0.d0
	do i=0,n-1
		dens=0.d0
		do j=-m,m
			dens=dens+f(i,j)                                  ! sum of f(x,v) over V
		enddo
		rho(i)=dens*dv                           ! rho: density at the X grid points
		densm=densm+dens
	enddo
	densm=densm*dv/dfloat(n)           ! mean density, or the contribution from ions
	do i=0,n-1
		rho(i)=rho(i)-densm    ! minus ions part, to ensure the mean density be zero
		d(i)=rho(i)*dx**2
	enddo

	! *************************************

	! Calculation for potential, use 'tridag' method
	! phi(i+1)-2*phi(i)+phi(i-1)=rho*dx**2
	! *************************************
	aa=1.0d0
	bb=-2.0d0
	cc=1.0d0
	bet=bb
	phi(0)=d(0)/bet
	do i=1,n-1
		! decomposition
		gam(i)=cc/bet
		bet=bb-aa*gam(i)
		! forward substitution
		phi(i)=(d(i)-aa*phi(i-1))/bet
	enddo
	do i=n-2,0,-1
		! backward substitution
		phi(i)=phi(i)-gam(i+1)*phi(i+1)
	enddo
	! *************************************

	! Calculation for EField, still use 'tridag' method
	! ef(i+1)+4*ef(i)+ef(i-1)=3*(phi(i+1)-phi(i-1))/dx
	! *************************************
	tddx=3.0d0/dx
	do i=1,n-2                           ! calculation for the right hand side (rhs)
		d(i)=(phi(i+1)-phi(i-1))*tddx
	enddo
	d(0)=(phi(1)-phi(n-1))*tddx
	d(n-1)=(phi(0)-phi(n-2))*tddx

	aa=1.0d0
	bb=4.0d0
	cc=1.0d0
	bet=bb
	ef(0)=d(0)/bet
	do i=1,n-1
		! Decomposition
		gam(i)=cc/bet
		bet=bb-aa*gam(i)
		! forward substitution
		ef(i)=(d(i)-aa*ef(i-1))/bet
	enddo
	do i=n-2,0,-1
		! backward substitution
		ef(i)=ef(i)-gam(i+1)*ef(i+1)
	enddo

	! Ensure Efield periodicity on boundaries and the mean field value equal to zero
	call efieldcorr
	! *************************************

end subroutine efield3

!*---------------------------------------------------------------------------------*
subroutine outputpara ! output parameters
	use variable
	open(unit=14,file='parameter.out',form='formatted',status='replace')
	write(14,*) 'n, m, nt, nplot, ndiag, dt, dx, dv, vmax, rk0, xl'
	write(14,'(5(I10/),6(F10.4/))') n, m, nt, nplot, ndiag, dt, dx, dv, vmax, rk0 &
	                                , xl
	close(14)
end subroutine outputpara

!*---------------------------------------------------------------------------------*
subroutine history                                             ! Diagnostics routine
	! We calculates:
	!           mean density: <N>
	!           mean momentum: <P>
	!           mean eneygies: 1. kinetic <Ekin>; 2. electric <Eelect>; 
	!                          3. total <Etotal>=<Ekin>+<Eelect>.
	!                                                                 L
	!                    where <g>=1/L*{sum over X of a function g(X)}
	!                                                                 0
	use variable
	implicit real*8 (A-H,O-Z)

	real*8,dimension(-1:n) :: ek, p


	! calculation for output
	! *************************************
	fn=0.d0
	eek=0.d0
	eef=0.d0
	pp=0.d0
	s1=0.d0
	s2=0.d0
	do i=0,n-1
		hh=0.d0
		h=0.d0
		do j=-m,m
			h=h+v(j)*f(i,j)                                   ! momentum calculation
			hh=hh+v(j)*v(j)*f(i,j)                      ! kinetic energy calculation
			fn=fn+f(i,j)                                             ! total density
			ff=abs(f(i,j))+1.0d-16
			s1=s1+ff*ff                                    ! entropy, quadratic form
!			s2=s2+ff*log(ff)                            ! entropy, logarithmic form
		enddo
		ek(i)=hh*dv                                                 ! kinetic energy
		p(i)=h*dv                                                         ! momentum
		eek=eek+ek(i)                              ! mean kinetic energy calculation
		pp=pp+p(i)                                       ! mean momentum calculation
		eef=eef+ef(i)*ef(i)                       ! mean electric energy calculation
	enddo
	fn=fn*dv/dfloat(n)                                                ! mean density
	pp=pp/dfloat(n)                                                  ! mean momentum
	eek=0.5d0*eek/dfloat(n)                                    ! mean kinetic energy
	eef=0.5d0*eef/dfloat(n)                                   ! mean electric energy
	s1=-s1*dx*dv                          ! s1 test value of entropy, quadratic form
	s2=-s2*dx*dv                        ! s2 test value of entropy, logatithmic form
	! *************************************

	if(it==1) et0=eek+eef                             ! Get the initial total energy

	! Output to an external file
	!   t: time; fn: mean density; pp: mean momentum; eek: kinetic energy; eef: mean
	!   electric energy; eek+eef: total energy; s1,s2: entropy
	write(15,'(1x,F7.3,7ES18.8)') t, fn, pp, eek, eef, eek+eef, s1, s2
	write(stdout,'(A3,I5,A3,F8.3,A15,ES15.8,A20,F7.3)') 'it=', it, 't=', t, &
			'Total energy=', eek+eef, 'Etotal change(%)=', (eek+eef-et0)/et0*100
	
end subroutine history

!*---------------------------------------------------------------------------------*
subroutine plotf                     ! A routine used to write F to an external file

	use variable
	implicit real*8 (A-H,O-Z)

	character*4,external :: num2str
	character*8 output
	output='fxv'//num2str(nint(t/dt))
	open(unit=16,file=output,form='formatted')

	do i=0,n-1
		do j=-m,m
			! output i, j, x, v, f(x,v)
			write(16,'(2I6,2F10.4,ES18.8)') i, j, i*dx, j*dv, f(i,j)
		enddo
	enddo

	close(16)
	
end subroutine plotf

!*---------------------------------------------------------------------------------*
subroutine plotef              ! A routine used to write E Field to an external file

	use variable
	implicit real*8 (A-H,O-Z)

	character*4,external :: num2str
	character*8 output
	output='efx'//num2str(nint(t/dt))
	open(unit=17,file=output,form='formatted')

	do i=0,n-1
		! output i, x, ef(x), rho(x)
		write(17,'(I6,F10.4,2ES18.8)') i, i*dx, ef(i), rho(i)
	enddo

	close(17)
	
end subroutine plotef

!*---------------------------------------------------------------------------------*
! Function num2str from jzhu pic1d.f90, modified by hsxie, 2011-08-26
function num2str(num) 
	! used to create filename for subroutine 'plotf' & 'plotef'
	character*4 num2str
	integer num
	integer,dimension(4) :: nn
	do i=1,4
		nn(i)=num/10**(4-i)
		num=mod(num,10**(4-i))
		nn(i)=nn(i)+48
		num2str(i:i)=char(nn(i))
	end do
	return
end

!***********************************************************************************
!*                     End of the Program & References                             *
!*                Part of this code is inherited from [1-3]                        *
!   [1] Fijalkow, E., A numerical solution to the Vlasov equation, Computer Physics 
! Communications, 1999, 116, 319 - 328.
!   [2] http://cpc.cs.qub.ac.uk/summaries/ADJQ_v1_0.html.
!   [3] PICES1D code by Zhihong Lin, University of California, Irvine, 2006, 
! http://phoenix.ps.uci.edu/zlin/pic1d/.
!***********************************************************************************