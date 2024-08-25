% 2013-03-03 17:03
% seems OK, not perfect
close all; clear; clc;
nt=3000; dt=5e-2; t0=20.0;
% tmp=2;
% nt=nt*tmp;dt=dt/tmp;

L=100; nz=500; dz=L/nz; z=0:dz:L;
Ex=zeros(1,nz+1); By=zeros(1,nz+1); ux=zeros(1,nz+1); 
za=6*L/10; zb=7*L/10; zc=8*L/10; zd=9*L/10;
n00=0.8; % omega_p=sqrt(n00)
nk=10.0;
n0=zeros(1,nz+1)+n00.*((z>za)&(z<zb))+...
    n00.*(nk*(z-zc)/(zd-zc)+1).*(abs(z-(zc+zd)/2)<(zd-zc)/2);
eom=1e0;

Amp=1; w=1.0;

h=figure('unit','normalized','Position',[0.01 0.47 0.6 0.45]);
set(gcf,'DefaultAxesFontSize',15);
subplot(331); plot(z,n0,'LineWidth',2); xlabel('z');ylabel('n_0(z)');
axis tight; grid on;
for it=1:nt
    
    for iz=2:nz+1
        Ex(iz)=Ex(iz)-dt*(By(iz)-By(iz-1))/dz+n0(iz)*ux(iz)*dt;
    end

    if(abs(it*dt-t0)<=t0)
        % Source Pulse
        pulse=Amp*sin((it*dt-t0)*w).*exp(-((it*dt-t0)/3)^2);
%         Ex(1)=Ex(1)+pulse;
        Ex(1)=pulse;
    else
        Ex(1)=Ex(2);
    end
%     Ex(nz+1)=Ex(nz);

    ux=ux-dt*Ex*eom; % cal ux
    
    for iz=1:nz
        By(iz)=By(iz)-dt*(Ex(iz+1)-Ex(iz))/dz;
    end
    
%     By(1)=By(2);
    By(nz+1)=By(nz);
    
    ymax=1.2*Amp;
    
    subplot(339);
    plot(z,Ex,'b',z,By,'g--',z,ux,'r-.','LineWidth',1.5);
    xlim([0,L]); ylim([-ymax,ymax]);
    rectangle('Position',[za, -ymax, zb-za, 2*ymax],'LineWidth',1.5);
    rectangle('Position',[zc, -ymax, zd-zc, 2*ymax],'LineWidth',1.5);
    text(0,ymax/2,['T=',num2str(it-1)]);
    
    
    nttmp=floor(nt/8);
    if(mod(it-1,nttmp)==0)
        subplot(3,3,floor((it-1)/nttmp)+2);
        plot(z,Ex,'b',z,By,'g--',z,ux,'r-.','LineWidth',1.5);
        xlim([0,L]); ylim([-ymax,ymax]);
        rectangle('Position',[za, -ymax, zb-za, 2*ymax],'LineWidth',1.5);
        rectangle('Position',[zc, -ymax, zd-zc, 2*ymax],'LineWidth',1.5);
        text(0,ymax/2,['T=',num2str(it-1)]);
    end
    
    pause(0.002);
end

str=['nt=',num2str(nt),',dt=',num2str(dt),',nz=',num2str(nz),',dz=',...
    num2str(dz),',n00=',num2str(n00),',eom=',num2str(eom)];
% title(str); 
subplot(332);legend('Ex','By','ux',3); legend('boxoff');

print('-dpng',[str,'.png']);

