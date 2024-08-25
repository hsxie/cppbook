% stellarator_shape_3d.m
% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2012-01-05 15:09
% Ref:
% Hirshman, S. & Lee, D. MOMCON: A spectral code for obtaining
% three-dimensional magnetohydrodynamic equilibria Computer Physics
% Communications, 1986, 39, 161 - 172 
clear;clc;

[theta,phi]=meshgrid(0:2*pi/50:2*pi);
zeta=6.*phi;

% From Eq(14)
Rb=1.72+...
   (0.214.*cos(theta)-0.0564.*cos(theta-zeta)+0.0035.*cos(theta+zeta))+...
   (0.00243.*cos(2.*theta)-0.00795.*cos(2.*theta-zeta));
Zb=(0.251.*sin(theta)+0.0634.*sin(theta-zeta)+0.0035.*sin(theta+zeta))+...
   (-0.00561.*sin(2.*theta)+0.00805.*sin(2.*theta-zeta));
X=Rb.*cos(phi);
Y=Rb.*sin(phi);
Z=Zb;

subplot(221);plot(Rb(1,:),Zb(1,:),'linewidth',2);xlabel('Rb');ylabel('Zb');title('At \phi=0');
subplot(222);mesh(X,Y,Z);xlabel('X');ylabel('Y');zlabel('Z');
axis equal;hidden off;colormap([0 1 1]);title('Stellatator grids'); box on;
subplot(223);mesh(X,Y,Z);view(0,90);xlabel('X');ylabel('Y');zlabel('Z');
axis equal;hidden off;colormap([0 1 1]);title('Stellatator grids');
subplot(224);mesh(X,Y,Z);view(90,0);xlabel('X');ylabel('Y');zlabel('Z');
axis equal;hidden off;colormap([0 1 1]);title('Stellatator grids');

