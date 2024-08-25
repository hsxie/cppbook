% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2012-12-20 20:00
% ideal Ballooning mode, initial value solver
% close to shooting method for unstable root
close all; clear; clc;
s=0.4; a=0.8; % s, alpha

nt=2000; dt=0.05;
dx=0.2; xmax=50; x=-xmax:dx:xmax; nx=length(x);
phi=zeros(nx,nt+1);
phi(:,1)=exp(-x.^2); phi(:,2)=phi(:,1);

dx2=dx*dx; dt2=dt*dt;
p=2.0.*(s.*x-a.*sin(x)).*(s-a.*cos(x));
q=a.*(cos(x)+(s.*x-a.*sin(x)).*sin(x));
r=(1.0+(s.*x-a.*sin(x)).^2);

set(gcf,'DefaultAxesFontSize',15);
for it=2:nt
    for j=2:nx-1
        phi(j,it+1)=(phi(j+1,it)-2*phi(j,it)+phi(j-1,it))*dt2/dx2+...
            p(j)/r(j)*(phi(j+1,it)-phi(j-1,it))*dt2/(2*dx)+...
            q(j)/r(j)*phi(j,it)*dt2+...
            2*phi(j,it)-phi(j,it-1);
    end
    if(mod(it,100)==0)
        subplot(211);plot(x,phi(:,it),'b','LineWidth',2);
        title(['s=',num2str(s),', \alpha=',num2str(a),...
            ', t=',num2str(it*dt)]);
        xlabel('\theta'); ylabel('\phi'); % xlim([0,xmax]);
    end
    pause(0.01);
end
tt=0:dt:nt*dt;
phi0=phi(floor(end/2),:);
subplot(223);plot(tt,phi0,'LineWidth',2); title('\phi(0,t)');
subplot(224);plot(tt,log(phi0),'LineWidth',2);
gamma=(log(phi0(end))-log(phi0(floor(end/2))))/(tt(end)-tt(floor(end/2)));
title(['ln[\phi(0,t)], \gamma=',num2str(gamma)]);

print(gcf,'-dpng',['s=',num2str(s),',alpha=',num2str(a),...
    ',gamma=',num2str(gamma),'.png']);