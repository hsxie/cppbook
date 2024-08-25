function kdv_MOL
close all; clear all; clf;
% Space coordinates
dx = 0.1;  x = (-8+dx:dx:8)';
nx = length(x);
% Time steps
k = dx^3; nsteps = 2.0 /k;
% Initial condition
% u = onesoliton(x,16,0);
 u =-8*exp(-x.^2);
%  u=-6./cosh(x).^2;
% u=onesoliton(x,16,0)+onesoliton(x,4,0);
% u=onesoliton(x,16,4)+onesoliton(x,4,-4)+onesoliton(x,8,-3);
set(gcf,'doublebuffer','on');
for ii=1:nsteps
    % Runge-Kutta step
    k1=k*kdvequ(u,dx);
    k2=k*kdvequ(u+k1/2,dx);
    k3=k*kdvequ(u+k2/2,dx);
    k4=k*kdvequ(u+k3,dx);
    u=u+k1/6+k2/3+k3/3+k4/6;
    % Animate every 10th step
    if mod(ii,10)==0
        plot(x,-u,'LineWidth',2); axis([-8,8,-2,12]); drawnow; pause(0.1);
    end
end
end

function u=onesoliton(x,v,x0)
    u=-v/2./cosh(.5*sqrt(v)*(x-x0)).^2;
end

function dudt=kdvequ(u,dx)
    % KdV equation: dudt = 6*u*dudx - d^3u/dx^3
    u = [u(end-1:end);u;u(1:2)];
    dudt = 6*(u(3:end-2)).*(u(4:end-1)-u(2:end-3))/2/dx - ...
           (u(5:end)-2*u(4:end-1)+2*u(2:end-3)-u(1:end-4))/2/dx^3;
end
