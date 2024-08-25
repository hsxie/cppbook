program planetorbit
  x0=1; vx0=0; y0=0; vy0 = 1
  read (*,*) dt
  N = 30/dt

  xel0=x0; vxel0=vx0; yel0=y0; vyel0=vy0 ! el=Euler

  xlf0=x0; vxlf0=vx0; ylf0=y0; vylf0=vy0 ! lf=Leapfrog
  xlf1=xlf0+vxlf0*dt; ylf1=ylf0+vylf0*dt
  xhlf0=(xlf0+xlf1)/2; yhlf0=(ylf0+ylf1)/2
  t=0
  do i = 0, N
    t=t+dt
    xel1 = xel0 + vxel0*dt
    yel1 = yel0 + vyel0*dt
    rel = sqrt(xel0*xel0 + yel0*yel0)
    fxel = -xel0/rel**3
    fyel = -yel0/rel**3
    vxel1 = vxel0 + fxel*dt
    vyel1 = vyel0 + fyel*dt
	
    xhlf1 = xhlf0+vxlf0*dt; 
	yhlf1 = yhlf0 + vylf0*dt;
    rlf = sqrt(xhlf0*xhlf0 + yhlf0 *yhlf0 )
    fxlf = -xhlf1/rlf**3
    fylf = -yhlf1/rlf**3
    vxlf1 = vxlf0 + fxlf*dt
    vylf1 = vylf0 + fylf*dt
    ! if(mod(i,N/10).eq.2)
    write(*,*) t, xel0, yel0, -1/rel+(vxel0*vxel0+vyel0*vyel0)/2,&
          xhlf0, yhlf0, -1/rlf+(vxlf0*vxlf0+vylf0*vylf0)/2

    xel0=xel1; yel0=yel1; vxel0=vxel1; vyel0=vyel1
    xhlf0=xhlf1; yhlf0=yhlf1; vxlf0=vxlf1; vylf0=vylf1
  enddo
end program
