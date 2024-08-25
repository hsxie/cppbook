% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2017-04-05 09:08
% Solve the 1D Burgers equation

close all; clear; clc;

nu=0.02; nt=1001; dt=0.02; dx=0.05; x=-8:dx:8; L=max(x)-min(x);
u0=exp(-(x+3).^2); u=u0; nx=length(x)-1; utmp=u;
figure('unit','normalized','position',[0.1,0.1,0.4,0.5],...
    'DefaultAxesFontSize',12);
for it=1:nt
    utmp(2:nx)=u(2:nx)-u(2:nx).*(u(3:(nx+1))-u(1:(nx-1)))*dt/(2*dx)+...
        nu*(u(3:(nx+1))-2*u(2:nx)+u(1:(nx-1)))*dt/(dx*dx);
    u=utmp;
    if(mod(it,floor(nt/5))==1)
        if(it<=1)
            plot(x,u,'r:','LineWidth',2);hold on;
        else
            plot(x,u,'b','LineWidth',2);hold on;
        end
        [ym,idx]=max(u);
        text(x(idx),ym+0.05,['t=',num2str((it-1)*dt)]);
    end
end
title(['Burgers, \nu=',num2str(nu),', L=',num2str(L),', dx=',num2str(dx),', dt=',num2str(dt),...
    ', nt=',num2str(nt)]);
xlim([min(x),max(x)]); ylim([0,1.1]); xlabel('x');ylabel('u');
print(gcf,'-dpng',['burgers_nu=',num2str(nu),',L=',num2str(L),...
    ',dx=',num2str(dx),',dt=',num2str(dt),...
    ',nt=',num2str(nt),'.png']);

