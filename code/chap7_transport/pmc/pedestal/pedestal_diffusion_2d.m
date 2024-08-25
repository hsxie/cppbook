% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2012-10-06 11:43
% pedestal_diffusion_2d.m
%
% Solution of the heat equation u_t=D(u_xx+u_yy), for SOL (Tokamak)
% using expectation of the heat kernel and random walk methods (MC)
%
close all; clear all; clc;
N=4e6; M=200; D=1; r0=1.0; 
T=0.0036; % 2*pi*L/c_s
xmax=2; xmin=-xmax; ymax=xmax; ymin=-xmax;
x=linspace(xmin,xmax,M); y=linspace(ymin,ymax,M);
dx=x(2)-x(1); dy=y(2)-y(1);
[xx,yy]=meshgrid(x,y);
% ue=0.5*((x>-2.5)&(x<-0.5))+((x>0.5)&(x<2.5));

figure('Unit','normalized','position',[0.1 0.65 0.55 0.25]);
set(gcf,'DefaultAxesFontSize',15);

ue=1.0*(xx.^2+yy.^2<=r0^2);
for k=1:N
    xi(k)=2*r0*(rand-0.5);yi(k)=2*r0*(rand-0.5);
    while((xi(k)^2+yi(k)^2)>=r0^2)
        xi(k)=2*r0*(rand-0.5);yi(k)=2*r0*(rand-0.5);
    end
    touti(k)=0; % record the time when a particle is at outer region
end
xii=xi; yii=yi;
m1=1/(N*dx); m2=1/(N*dy);
u=fxy(xi,yi,xx,yy,dx,dy,N);
subplot(131);pcolor(xx,yy,u); shading interp;
title('contour for t=0');xlabel('R-R0');ylabel('Z');
subplot(133);
xmid=xx(round(end/2),:); umid0=u(round(end/2),:);
plot(xmid,umid0,':','LineWidth',2);
dt=0.002;
for t=1:50
    xi=xii+sqrt(2*D*dt*t)*randn(1,N);
    yi=yii+sqrt(2*D*dt*t)*randn(1,N);
    for k=1:N % outer region of SOL
        if((xi(k)^2+yi(k)^2)>=r0^2)
            touti(k)=touti(k)+dt;
        else
            touti(k)=0;
        end
        randt=T*rand;
        if(touti(k)>randt)
            xi(k)=NaN; yi(k)=NaN;
        end
    end
    u=fxy(xi,yi,xx,yy,dx,dy,N);
    subplot(132);pcolor(xx,yy,u); shading interp;
    title(['N=',num2str(N,'%3.1e'),', t=',num2str(t*dt)]);
    xlabel('R-R0');ylabel('Z');
    subplot(133);
    umid=u(round(end/2),:);
    plot(xmid,umid0,':',xmid,umid,'ro','LineWidth',2);
    xlabel('R-R0');ylabel('Density');
    title(['D=',num2str(D),', T=',num2str(T),' (Z=0)']);
    pause(0.1);
    drawnow;
end
