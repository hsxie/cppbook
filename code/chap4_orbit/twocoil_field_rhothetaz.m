close all;clear;clc;
c0=2e-7; % mu0/2pi=2e-7
I=1e6; % current, unit: A
a=0.2; L=2.0;
eps=1e-4;
[X,Y,Z]=meshgrid(-2*a:a/10:2*a,-2*a:a/10:2*a,-1.5*L:L/10:1.5*L);
RHO2=X.^2+Y.^2; RHO=sqrt(RHO2);
THETA=atan2(Y,X);

k2=4*a.*RHO./((a+RHO).^2+Z.^2);
k2=k2-eps*(k2==1);
[K,E]=ellipke(k2); % 0<=k2<=1
a1=sqrt((a+RHO).^2+Z.^2);
b1=(a-RHO).^2+Z.^2;
b1=b1+eps*(b1==0);
b2=a^2+RHO.^2+Z.^2;
b3=a^2-RHO.^2-Z.^2;
Brho=c0*I.*Z./((RHO+eps*(RHO==0)).*a1).*(b2./b1.*E-K);
Bz=c0*I.*Z./a1.*(b3./b1.*E+K);
Btheta=0.*Bz;

Bx=Brho.*cos(THETA)-Btheta.*sin(THETA);
By=Brho.*sin(THETA)+Btheta.*cos(THETA);

[rho,theta,z]=meshgrid(0.01*a:0.25*a:0.9*a,0:pi/6:2*pi,0.01);
sx=rho.*cos(theta); sy=rho.*sin(theta); sz=0.0.*sx+0.01;
streamline(X,Y,Z,Bx,By,Bz,sx,sy,sz);zlim([-L,L]);
xlabel('x');ylabel('y');zlabel('z');axis equal;
view(-8,38);
hold on;
B=sqrt(Bx.^2+By.^2+Bz.^2);
Zmin=min(Z(:));
% hcont=contourslice(X,Y,Z,B,[],[],-L:0.4*L:L);
hcont=contourslice(X,Y,Z,B,[],[],0.1*L);
set(hcont,'EdgeColor','red','LineWidth',.5);

