close all;clear;clc;
R0=3.0;a=R0/4;r0=0.8*a;dr=0.3*r0;
[theta,phi]=meshgrid(0:2*pi/80:2*pi,0*pi:2*pi/80:2.0*pi);

m=2;n=1;

r=r0+dr.*cos(m.*theta-n.*phi);
R=R0+r.*cos(theta);
Z=r.*sin(theta);
X=R.*cos(phi);
Y=R.*sin(phi);
figure; set(gcf,'DefaultAxesFontSize',15);
surf(X,Y,Z);
axis equal; 
% shading interp;
% colormap(jet); alpha(0.5);
str=['(m=',num2str(m),', n=',num2str(n),') kink mode'];
title(str);
print('-dpng',str);

%% cylinder
% close all;clear;clc;
% R0=3.0;a=R0/2;r0=0.5*a;dr=0.3*r0;
% [theta,phi]=meshgrid(0:2*pi/50:2*pi,0*pi:2*pi/50:2.0*pi);
% r=r0+dr.*cos(2.*theta-phi);
% X=r.*cos(theta);
% Y=r.*sin(theta);
% Z=R0.*phi;
% surf(X,Y,Z);
% axis equal;