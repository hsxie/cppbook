% draw_field_line.m
close all;clear;clc;
R0=3; a=R0*0.3; n=1; m=3; q=m/n;
[theta,phi]=meshgrid(0:2*pi/50:2*pi,0:2*pi/50:1.5*pi);
X=(R0+a.*cos(theta)).*cos(phi);
Y=(R0+a.*cos(theta)).*sin(phi);
Z=a.*sin(theta); 
surf(X,Y,Z); hold on;
hidden off;
colormap hsv;
shading interp;
alpha(0.2);
axis equal;
% axis('equal','square','off');

theta2=0:pi/1000:m*n*2*pi;
phi2=q.*theta2;
strtitle=['B Field Line, q=m/n=',num2str(m),'/',num2str(n),', R0=',...
    num2str(R0),'m, a=',num2str(a),'m'];
x=(R0+a.*cos(theta2)).*cos(phi2);
y=(R0+a.*cos(theta2)).*sin(phi2);
z=a.*sin(theta2);
plot3(x,y,z,'g','LineWidth',2);alpha(0.2);xlabel('x');ylabel('y');zlabel('z');
title(strtitle);
