close all;clear;clc;
R0=3.0;a=R0/4;r0=0.8*a;dr=0.1*r0;
[theta,phi]=meshgrid(0:2*pi/150:2*pi,-0.5*pi:2*pi/60:1.2*pi);

m=15;n=10;

envelope=exp(-0.8.*(abs(theta-pi)-pi).^2);

tmp=(theta-pi).^2-pi^2;

r=r0+(dr.*envelope).*cos(m.*theta-n.*phi);
R=R0+r.*cos(theta);
Z=r.*sin(theta);
X=R.*cos(phi);
Y=R.*sin(phi);
figure; set(gcf,'DefaultAxesFontSize',15);
surf(X,Y,Z);
axis equal; 
% shading interp;
% colormap(jet); alpha(0.5);
str=['(m=',num2str(m),', n=',num2str(n),') ballooning mode'];
title(str);
print('-dpng',str);

