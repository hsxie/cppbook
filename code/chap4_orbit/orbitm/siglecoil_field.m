clear;clc;
a=0.35;L=2;
[X,Y]=meshgrid(-1:0.01:1);
XY=X.^2+Y.^2;
Xa=abs(X); Ya=abs(Y);
k2=4*a*Xa./(a^2+XY+2*a*Xa);
[K,E]=ellipke(k2); % 0<=k2<=1
aXY=a^2+XY-2*a*Xa;
aXY=aXY+eps*(aXY==0);
aE1=(a^2-XY)./aXY.*E;
aE2=(a^2+XY)./aXY.*E;
By=1/pi./sqrt(a^2+XY+2*a*Xa).*(aE1+K);
Bx=1/pi*Y./(X+eps*(X==0))./sqrt(a^2+XY+2*a*Xa).*(aE2-K);
vy=0; vx=[1e-4:a/30:a];
[Vx, Vy]=meshgrid(vx,vy);
streamline(X,Y,Bx,By,Vx,Vy);
% hold on;
% streamline(-X,Y,-Bx,By,-Vx,Vy);
% streamline(-X,-Y,-Bx,-By,-Vx,-Vy);
% streamline(X,-Y,Bx,-By,Vx,-Vy);