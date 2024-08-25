% tokamak_shape_3d.m
% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2012-01-05 15:49
% Ref:
% White, R. B. The Theory of Toroidally Confined Plasmas World Scientific
% Imperial College Press, 2001, P143

clear;clc;
R0=3;a=R0/3;kappa=1.0;delta=0.0;b=0.0;

[theta,phi]=meshgrid(0:2*pi/50:2*pi,0:2*pi/50:1.5*pi);
X=(R0-b+(a+b.*cos(theta)).*cos(theta+delta.*sin(theta))).*cos(phi);
Y=(R0-b+(a+b.*cos(theta)).*cos(theta+delta.*sin(theta))).*sin(phi);
Z=(kappa*a).*sin(theta); 
mesh(X,Y,Z);
axis equal;
hidden off;
colormap([0 1 1]);
title(['Tokamak grids, with R0=',num2str(R0),'m ,a=',num2str(a),'m ,b=',...
    num2str(b),'m ,\kappa=',num2str(kappa),' ,\delta=',num2str(delta)]);