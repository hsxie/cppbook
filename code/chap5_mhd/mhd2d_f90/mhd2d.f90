!********************************************************************************
!*   mhd2d.f90, 2D Compressible MHD Code for Tearing Mode, 2012-06-17 11:55     *
!*                 Ref: [Fu1995], section 7.4 and Appendix 2                    *
!*                                                                              *
!*         Hua-sheng XIE, IFTS-ZJU, huashengxie@{gmail.com, zju.edu.cn}         *
!*                                                                              *
!*         Boundary condition: x -- ux=bx=0, p_x=uz_x=rho_x=0                   *
!*                             z -- periodic boundary                           *
!*                                                                              *
!*  [Fu1995] 傅竹风 & 胡友秋, 空间等离子体数值模拟, 安徽科学技术出版社, 1995    *
!********************************************************************************
!            
!           ↑ x
!           |
!         Lx|------------------------------  
!           |                              |
!           |                              |
!          0|                              |Lz
!     -------------------------------------------> z
!           |                              |
!           |                              |
!           |                              |
!        -Lx|------------------------------   
!           |
!            

Module variable

	implicit none
	
	! grid parameters
	integer, parameter :: nii=5, njj=5
	integer, parameter :: ni=2**nii+1                         ! ni, x grid number
	integer, parameter :: nj=2**njj                           ! nj, z grid number
	integer, parameter :: nt=5000
	real*8, parameter :: dx=1.0, dz=2.0
	real*8 :: dt=0.05     ! giving a dt < min(dx,dz)/[sqrt(1.0+0.5*gamma*beta)*va]
	real*8 :: rdx=0.5/dx, rdz=0.5/dz, rdx2=1.0/(dx*dx), rdz2=1.0/(dz*dz)

	! physical parameters
	real*8, parameter :: gamma=1.66667              ! parameters in mhd equations
	real*8, parameter :: eta=0.01
	real*8, parameter :: nu=0.05
	real*8, parameter :: beta=0.5                                  ! beta in x=Lx
	real*8, parameter :: rho0=1.0             ! rho0 -- mass density rho0 in x=Lx
	real*8, parameter :: b0=1.0                                        ! b0 -- B0
	real*8, parameter :: bl=3.0                    ! bl -- width of current sheet
	real*8 :: va                                                ! Alfven velocity
	real*8 :: p0                                      ! p0 -- pressure p0 in x=Lx
	real*8 :: t0                                   ! t0 -- temperature T0 in x=Lx
	real*8 :: bxm                                           ! bxm -- log(max(Bx))
	
	! other parameters
	integer, parameter :: nplot=(nt+1)/50, ndiag=1
	real*8, parameter :: pi=3.14159265358979, small=1.0e-10, large=1.0e10
	real*8, parameter :: cfl=1.5  ! Courant–Friedrichs–Lewy condition parameter
	integer, parameter ::  nw=1            ! nw -- number of initial perturbation
	real*8, dimension(10) :: am                            ! am -- prtb amplitude
	data am /0.01, 9*0.0/
	real*8, parameter :: abd=10.0                ! abd -- perturbation width in x
	integer ::  kw=0                     ! kw -- warning, =0 normal, =1 divergent

	real*8, dimension(ni,nj,5) :: x                ! 1 2 3 4 5 --> rho p ux uz ay
	real*8, dimension(ni,nj,5) :: y, f1, f2, f3, f4         ! temp arrays for rk4
	real*8, dimension(ni,nj) :: bx, bz, pt                           ! pt=p+b^2/2
	
	integer :: stdout=1      !stdout=0: output to display; 1: to file "mhd2d.out"
	
end Module variable

!********************************************************************************
Program mhd2d

	use variable
	implicit none
	integer cause,it,i,j
	real*8 t
	real ctime0, ctime1

	call CPU_TIME(ctime0)                              ! obtain starting CPU time
	if(stdout==1) open(stdout,file='mhd2d.out',status='replace')
	write(stdout,*)'ntime       log(max(Bx))'

	open(15,file='bxm.dat',status='replace')
	
	call initial
	do it=1, nt

		if(mod(it-1,nplot)==0) then
			call plot(it)
		endif
		
		call bxmax
		write(15,'(I6,F10.4,ES18.8)') it, it*dt, bxm
		if(mod(it-1,ndiag)==0) then
			write(stdout,*) it, log(bxm)
		endif

		t=t+dt
		call rk4
		do i=1,ni
			do j=1,nj
				cause=1
				if(x(i,j,1)<0)	exit
				cause=2
				if(x(i,j,2)<0)	exit
			enddo
		enddo
		cause=3
		if(kw==1)	exit
		cause=0
		
	enddo
	call exitinfo(cause,it)

	call CPU_TIME(ctime1)                                ! obtain ending CPU time
	write(stdout,*) 'Program CPU TIME=', ctime1-ctime0

	close(15)
	
end Program mhd2d

!********************************************************************************
Subroutine initial

	use variable
	implicit none
	integer i,j,njh,m
	real*8 s,s1,b,p,rho,sii,sij
	
	va=dsqrt(b0*b0/rho0)
	p0=0.5*beta*b0*b0
	t0=0.5*beta*va*va
	
	x=0.0
	do i=1,ni
		s=(i-ni)*dx/bl
		s1=b0*bl*log(cosh(s))
		b=b0*tanh(s)
		p=p0+0.5*(b0**2-b**2)
		rho=p/t0
		do j=1,nj
			x(i,j,5)=s1
			x(i,j,1)=rho
			x(i,j,2)=p
		enddo
	enddo
	njh=0.5*nj
	do m=1,nw
		do i=2,ni
			s=(i-ni)*dx/abd
			sii=exp(-s*s)*am(m)*b0
			do j=1,nj
				sij=sin((2.0*m*(j-njh)/(nj-1.0)+0.5)*pi)*dz*(nj-1)/(2.0*m*pi)
				x(i,j,5)=x(i,j,5)+sij*sii
			enddo
		enddo
	enddo
end Subroutine initial

!********************************************************************************
Subroutine rk4	                              ! 4-th Runge-Kutta time integration
	! required subroutine right(xo,xi)
	!      xi, variable; xo, right side of mhd equations
	! input: x -- variable array
	! output: x -- after one time step
	!         kw -- warning, =0 normal, =1 divergent
	
	use variable
	implicit none
	integer i,j,k
	
	call right(f1,x)	
	y=x+0.5*dt*f1
	call right(f2,y)
	y=x+0.5*dt*f2
	call right(f3,y)
	y=x+dt*f3
	call right(f4,y)
	x=x+dt*(f1+2.0*f2+2.0*f3+f4)/6.0
	
	do k=1,5                                             ! judge divergent or not
		do j=1,nj
			do i=1,ni
				if(abs(x(i,j,k))>large) then
					kw=1
					return
				endif
			enddo
		enddo
	enddo
	
	return
end Subroutine rk4

!********************************************************************************
Subroutine right(xo,xi)               ! for rk4, right hand side of mhd equations

	use variable
	implicit none
	integer i,j,k,jp1,jm1
	real*8, dimension(ni,nj,5) :: xi,xo            ! 1 2 3 4 5 --> rho p ux uz ay
	real*8, dimension(ni,nj) :: rrho                                      ! 1/rho

	! calculate Bx, Bz
	call calcbxz(xi)

	do i=1,ni
		do j=1,nj
			rrho(i,j)=1.0/xi(i,j,1)                                 ! rrho=1/rho
			pt(i,j)=xi(i,j,2)+0.5*(bx(i,j)**2+bz(i,j)**2)           ! pt=p+b^2/2
		enddo
	enddo

	do j=1,nj
		jp1=j+1
		jm1=j-1
		if(j==nj) jp1=1
		if(j==1) jm1=nj
		do i=2,ni-1
			xo(i,j,1)=-xi(i,j,3)*rdx*(xi(i+1,j,1)-xi(i-1,j,1))&
				-xi(i,j,4)*rdz*(xi(i,jp1,1)-xi(i,jm1,1))&
				-xi(i,j,1)*(rdx*(xi(i+1,j,3)-xi(i-1,j,3))&
					+rdz*(xi(i,jp1,4)-xi(i,jm1,4)))
			xo(i,j,2)=-xi(i,j,3)*rdx*(xi(i+1,j,2)-xi(i-1,j,2))&
				-xi(i,j,4)*rdz*(xi(i,jp1,2)-xi(i,jm1,2))&
				-gamma*xi(i,j,2)*(rdx*(xi(i+1,j,3)-xi(i-1,j,3))&
					+rdz*(xi(i,jp1,4)-xi(i,jm1,4)))
			xo(i,j,3)=-xi(i,j,3)*rdx*(xi(i+1,j,3)-xi(i-1,j,3))&
				-xi(i,j,4)*rdz*(xi(i,jp1,3)-xi(i,jm1,3))&
				+rrho(i,j)*((bx(i,j)*rdx*(bx(i+1,j)-bx(i-1,j))&
					+bz(i,j)*rdz*(bx(i,jp1)-bx(i,jm1))&
							-rdz*(pt(i+1,j)-pt(i-1,j)))&
				+nu*(rdx2*(x(i+1,j,3)+xi(i-1,j,3)-2.0*xi(i,j,3))&
					+rdz2*(xi(i,jp1,3)+xi(i,jm1,3)-2.0*xi(i,j,3))))
			xo(i,j,4)=-xi(i,j,3)*rdx*(xi(i+1,j,4)-xi(i-1,j,4))&
				-xi(i,j,4)*rdz*(xi(i,jp1,4)-xi(i,jm1,4))&
				+rrho(i,j)*((bx(i,j)*rdx*(bz(i+1,j)-bz(i-1,j))&
					+bz(i,j)*rdz*(bz(i,jp1)-bz(i,jm1))&
							-rdz*(pt(i,jp1)-pt(i,jm1)))&
				+nu*(rdx2*(xi(i+1,j,4)+xi(i-1,j,4)-2.0*xi(i,j,4))&
					+rdz2*(xi(i,jp1,4)+xi(i,jm1,4)-2.0*xi(i,j,4))))
			xo(i,j,5)=-xi(i,j,3)*rdx*(xi(i+1,j,5)-xi(i-1,j,5))&
				-xi(i,j,4)*rdz*(xi(i,jp1,5)-xi(i,jm1,5))&
				+eta*(rdx2*(xi(i+1,j,5)+xi(i-1,j,5)-2.0*xi(i,j,5))&
					+rdz2*(xi(i,jp1,5)+xi(i,jm1,5)-2.0*xi(i,j,5)))
		enddo
	enddo
	
	do j=1,nj                          ! z direction, periodic boundary condition
		jp1=j+1
		jm1=j-1
		if(j==nj) jp1=1
		if(j==1) jm1=nj
		xo(ni,j,1)=-xi(ni,j,4)*rdz*(xi(ni,jp1,1)-xi(ni,jm1,1))&
			-xi(ni,j,1)*(rdx*2*(-xi(ni-1,j,3))&
				+rdz*(xi(ni,jp1,4)-xi(ni,jm1,4)))
		xo(ni,j,2)=-xi(ni,j,4)*rdz*(xi(ni,jp1,2)-xi(ni,jm1,2))&
			-gamma*xi(ni,j,2)*(rdx*2*(-xi(ni-1,j,3))&
				+rdz*(xi(ni,jp1,4)-xi(ni,jm1,4)))
		xo(ni,j,3)=0.0
		xo(ni,j,4)=-xi(ni,j,4)*rdz*(xi(ni,jp1,4)-xi(ni,jm1,4))&
			+rrho(ni,j)*((bx(ni,j)*rdx*2*(-bz(ni-1,j))&
						-rdz*(pt(ni,jp1)-pt(ni,jm1)))&
			+nu*(rdx2*2*(xi(ni-1,j,4)-xi(ni,j,4))&
				+rdz2*(xi(ni,jp1,4)+xi(ni,jm1,4)-2.0*xi(ni,j,4))))
		xo(ni,j,5)=-xi(ni,j,4)*rdz*(xi(ni,jp1,5)-xi(ni,jm1,5))&
			+eta*(rdx2*2*(xi(ni-1,j,5)-xi(ni,j,5))&
				+rdz2*(xi(ni,jp1,5)+xi(ni,jm1,5)-2.0*xi(ni,j,5)))
	enddo

	do j=1,nj         ! x direction boundary condition, ux=bx=0, p_x=uz_x=rho_x=0
		xo(1,j,3)=0
		xo(1,j,5)=0
		
		xo(1,j,1)=xo(2,j,1)
		xo(1,j,2)=xo(2,j,2)
		xo(1,j,4)=xo(2,j,4)
	enddo
	
	return
end Subroutine right

!********************************************************************************
Subroutine calcbxz(xi)                          ! calculate Bx, Bz, need check!!!

	use variable
	implicit none
	integer i,j,k,jp1,jm1
	real*8, dimension(ni,nj,5) :: xi               ! 1 2 3 4 5 --> rho p ux uz ay
	
	! calculate Bx, Bz, need check!!!
	do i=2,ni-1
		do j=1,nj
			jp1=j+1
			if(j==nj) jp1=1
			jm1=j-1
			if(j==1) jm1=nj
			bx(i,j)=-rdz*(xi(i,jp1,5)-xi(i,jm1,5))       ! rdx=0.5/dx, rdz=0.5/dz
			bz(i,j)= rdx*(xi(i+1,j,5)-xi(i-1,j,5))
		enddo
	enddo
	! upper boundary, i=1
	do j=1,nj
		jp1=j+1
		if(j==nj) jp1=1
		jm1=j-1
		if(j==1) jm1=nj
		bx(1,j)=-rdz*(xi(1,jp1,5)-xi(1,jm1,5))
		bz(1,j)=2.0*rdx*(xi(2,j,5)-xi(1,j,5))
	enddo
	! low boundary, i=ni, symmetric
	do j=1,nj
		jp1=j+1
		if(j==nj) jp1=1
		jm1=j-1
		if(j==1) jm1=nj
		bx(ni,j)=-rdz*(xi(ni,jp1,5)-xi(ni,jm1,5))
		bz(ni,j)=0.0
	enddo

end Subroutine calcbxz

!********************************************************************************
Subroutine bxmax                                               ! calculate Bx max

	use variable
	implicit none
	integer i,j
	real*8 sss
	
	call calcbxz(x)                                                ! Calculate Bx
	
	bxm=-100.0
	do i=1,ni
		do j=1,nj
			sss=abs(bx(i,j))
			if(sss>bxm) bxm=sss
		enddo
	enddo
	if(bxm<small) bxm=small

end Subroutine bxmax

!********************************************************************************
Subroutine exitinfo(cause,it)

	use variable
	implicit none

	integer cause, it
	select case(cause)
		case(1)
			write(stdout,*) 'Negative density at',it,'-th step'
		case(2)
			write(stdout,*) 'Negative perssure at',it,'-th step'
		case(3)
			write(stdout,*) 'Divergent at',it,'-th step'
		case default
			write(stdout,*) 'it = nt_end=',it,'-th step'
	end select
end Subroutine exitinfo

!********************************************************************************
subroutine plot(it)         ! Routine used to write contour data to external file

	use variable
	integer i,j,it

	character*6,external :: num2str
	character*9 output1,output2,output3
	output1='rho'//num2str(it-1)
	output2='ux'//num2str(it-1)
	output3='Bx'//num2str(it-1)
	open(unit=16,file=output1,form='formatted')
	open(unit=17,file=output2,form='formatted')
	open(unit=18,file=output3,form='formatted')

	call calcbxz(x)                                                ! Calculate Bx
	
	do i=1,ni
		do j=1,nj
			! output i, j, x, z, -- rho(x,z), ux(x,z), Bx(x,z)
			write(16,'(2I6,2F10.4,ES18.8)') i, j, (i-(ni+1)/2)*dx, j*dz, x(i,j,1)
			write(17,'(2I6,2F10.4,ES18.8)') i, j, (i-(ni+1)/2)*dx, j*dz, x(i,j,3)
			write(18,'(2I6,2F10.4,ES18.8)') i, j, (i-(ni+1)/2)*dx, j*dz, bx(i,j)
		enddo
	enddo

	close(16)
	close(17)
	close(18)
	
end subroutine plot

!********************************************************************************
function num2str(num) 
	! Function num2str from jzhu pic1d.f90, modified by hsxie, 2011-08-26
	! used to create filename for subroutine 'plot*'
	character*6 num2str
	integer num
	integer,dimension(6) :: nn
	do i=1,6
		nn(i)=num/10**(6-i)
		num=mod(num,10**(6-i))
		nn(i)=nn(i)+48
		num2str(i:i)=char(nn(i))
	end do
	return
end