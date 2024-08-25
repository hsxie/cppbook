% tokamak_shape_3d.m
% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2012-08-31 10:52
% Ref:
% White, R. B. The Theory of Toroidally Confined Plasmas World Scientific
% Imperial College Press, 2001, P143

close all;clear;clc;
R0=3;a=0.9;kappa=1.5;delta=0.5;b=0.0;

[theta,phi]=meshgrid(0:2*pi/500:2*pi,-0.4*pi:2*pi/500:1.0*pi);
X=(R0-b+(a+b.*cos(theta)).*cos(theta+delta.*sin(theta))).*cos(phi);
Y=(R0-b+(a+b.*cos(theta)).*cos(theta+delta.*sin(theta))).*sin(phi);
Z=(kappa*a).*sin(theta);
R=R0+sqrt(X.^2+Y.^2);

h=surf(X,Y,Z,R); alpha(0.5); shading interp;
colormap(jet); 
% whitebg('b');

set(h,'FaceLighting','flat','FaceColor','interp',...
      'AmbientStrength',0.3);
camlight('right');
lightangle(45,60);
light('position',[-3,-1,3],'style','local');
material shiny;

axis equal; axis tight; axis off; hidden off;

n=1; m=3; q=m/n;
strtitle=['B Field Line, q=m/n=',num2str(m),'/',num2str(n),', R0=',...
    num2str(R0),'m, a=',num2str(a),'m, \kappa=',num2str(kappa),...
    ', \delta=',num2str(delta)];
theta2=0:2*pi/200:2*max(m,n)*pi;
phi2=q.*theta2;
x=(R0-b+(a+b.*cos(theta2)).*cos(theta2+delta.*sin(theta2))).*cos(phi2);
y=(R0-b+(a+b.*cos(theta2)).*cos(theta2+delta.*sin(theta2))).*sin(phi2);
z=(kappa*a).*sin(theta2);
hold on; plot3(x,y,z,'r','LineWidth',1);

title(strtitle);
print('-dpng',['B Field Line, m',num2str(m),'n=',num2str(n),', R0=',...
    num2str(R0),'m, a=',num2str(a),'m, kappa=',num2str(kappa),...
    ', delta=',num2str(delta),'.png']);