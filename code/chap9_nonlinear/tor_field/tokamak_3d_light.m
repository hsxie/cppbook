close all; clear; clc;
R0=3;
data=[R0/3,R0/5,R0/10; % a
    1.5, 1.5, 1.5;
    0.5, 0.5, 0.5;
    0.0, 0.0, 0.0];

nflux=size(data,2); h=figure; for iflux=1:nflux

    % a=R0/3;kappa=1.5;delta=0.5;b=0.0;
    a=data(1,iflux);kappa=data(2,iflux);
    delta=data(3,iflux);b=data(4,iflux);

    [theta,phi]=meshgrid(0:2*pi/50:2*pi,(-0.1-iflux*0.2)*pi:2*pi/50:(0.7+iflux*0.2)*pi);
    X=(R0-b+(a+b.*cos(theta)).*cos(theta+delta.*sin(theta))).*cos(phi);
    Y=(R0-b+(a+b.*cos(theta)).*cos(theta+delta.*sin(theta))).*sin(phi);
    Z=(kappa*a).*sin(theta);
    R=R0+sqrt(X.^2+Y.^2);

    surf(X,Y,Z,R); hold on;

end

camlight('right'); lightangle(45,60);
light('position',[-3,-1,3],'style','local'); material shiny;

axis equal; axis tight; axis off; hidden off;

strtitle=['Tokamak flux surfaces, R0=',...
    num2str(R0),'m, \kappa=',num2str(kappa),...
    ', \delta=',num2str(delta)];
title(strtitle);