close all;clear;clc;
R0=3.0;a=R0/3;r0=0.8*a; n=1; m=3;
[theta1,theta2]=meshgrid(-0.15*pi:2*pi/200:1.75*pi,0*pi:2*pi/20:2.0*pi);
q=m/n;
phi1=q.*theta1;
x0=(R0+r0.*cos(theta1)).*cos(phi1);
y0=(R0+r0.*cos(theta1)).*sin(phi1);
z0=r0.*sin(theta1);
RR0=sqrt(x0.^2+y0.^2);

aa=0.8*r0;bb=0.4*aa;
R=aa.*cos(theta2).*sin(phi1)-bb.*sin(theta2).*cos(phi1)+RR0;
Z=aa.*cos(theta2).*cos(phi1)+bb.*sin(theta2).*sin(phi1)+z0;
X=R.*cos(phi1);
Y=R.*sin(phi1);

figure; set(gcf,'DefaultAxesFontSize',15);
surf(X,Y,Z);
axis equal;
% shading interp;
% colormap(jet); alpha(0.5);
str=['(m=',num2str(m),', n=',num2str(n),') magnetic island'];
title(str);
print('-dpng',str);
