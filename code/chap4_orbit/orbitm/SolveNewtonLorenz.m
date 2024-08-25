function ydot=SolveNewtonLorenz(t,y)
global c q m Bfield rel Efield;

vsq = y(4)^2 + y(5)^2 + y(6)^2;

if (rel==1) % relativistic factor
    gamma = 1.0/sqrt(1 - vsq/c^2);
else % non-relativistic
    gamma=1;
end

Bvec = Bfield(y(1),y(2),y(3),t);
Bx=Bvec(1);By=Bvec(2);Bz=Bvec(3);
Evec=Efield(y(1),y(2),y(3),t);
Ex=Evec(1);Ey=Evec(2);Ez=Evec(3);

fac=q/(m*gamma);

ydot=[y(4);y(5);y(6);
   fac*(Ex+Bz*y(5)-By*y(6));
   fac*(Ey+Bx*y(6)-Bz*y(4));
   fac*(Ez+By*y(4)-Bx*y(5))];
