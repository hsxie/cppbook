close all;clear all;clc;
c0=1e-7; % mu0/4pi=1e-7
I=1e6; % current, unit: A
a=0.4; L=2.0;
eps=1e-4;
[RHO,Z]=meshgrid(0:0.01:1,-5:0.01:5);
k2=4*a.*RHO./((a+RHO).^2+Z.^2);
k2=k2-eps*(k2==1);
[K,E]=ellipke(k2); % 0<=k2<=1
a1=sqrt((a+RHO).^2+Z.^2);
b1=(a-RHO).^2+Z.^2;
b1=b1+eps*(b1==0);
b2=a^2+RHO.^2+Z.^2;
b3=a^2-RHO.^2-Z.^2;
Brho=c0*I.*Z./((RHO+eps*(RHO==0)).*a1).*(K-b2./b1.*E);
Bz=c0*I.*Z./a1.*(K+b3./b1.*E);

subplot(2,2,1);surf(RHO,Z,Brho);
subplot(2,2,2);surf(RHO,Z,Bz);
subplot(2,2,3);streamslice(RHO,Z,Brho,Bz);axis tight;
% subplot(2,2,4);streamline(RHO,Z,Brho,Bz);

% [X,Z]=pol2cart(RHO,Z);
% [Bx,Bz]=pol2cart(Brho,Bz);

% subplot(2,2,4);streamslice(X,Z,Bx,Bz);


vrho=[1e-4:a/30:a]; vz=0;
[Vrho, Vz]=meshgrid(vrho,vz);
subplot(2,2,4);streamline(RHO,Z,Brho,Bz,Vrho,Vz);

% [X,Y]=meshgrid(-1:0.01:1);
% XY=X.^2+Y.^2;
% Xa=abs(X); Ya=abs(Y);
% k2=4*a*Xa./(a^2+XY+2*a*Xa);
% [K,E]=ellipke(k2); % 0<=k2<=1
% aXY=a^2+XY-2*a*Xa;
% aXY=aXY+eps*(aXY==0);
% aE1=(a^2-XY)./aXY.*E;
% aE2=(a^2+XY)./aXY.*E;
% By=1/pi./sqrt(a^2+XY+2*a*Xa).*(aE1+K);
% Bx=1/pi*Y./(X+eps*(X==0))./sqrt(a^2+XY+2*a*Xa).*(aE2-K);
% vy=0; vx=[1e-4:a/30:a];
% [Vx, Vy]=meshgrid(vx,vy);
% streamline(X,Y,Bx,By,Vx,Vy);
