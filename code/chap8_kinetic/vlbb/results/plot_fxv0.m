% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2011-06-23 05:35
% plot_fxv.m, for data analysis of vlasov es 1d code vl1dper.f90
close all;clear;clc;

% modify here base on parameters of vl1dper.f90
vmax=8.0;rk0=.3;xl=2.0*pi/rk0; dt=0.125;
M=127;N=32;nplot=40;nt=1200;

dx=xl/(N-1);dv=vmax/M;
vv=-vmax:dv:vmax;
[ax,av]=meshgrid(0:dx:xl,-vmax:dv:vmax);
figure('position',[100 100 800 600]);
set(gca,'nextplot','replacechildren');
j=1;
for i=0:nplot:nt
    data=load(['fxv',num2str(i,'%4.4d')]);
    x=data(:,1);
    v=data(:,2);
    f=reshape(data(:,3),2*M+1,[]);
    if(i==0) fv0=mean(f'); maxfv0=max(fv0);
    end
    fv=mean(f');
    color=data(:,3);
%     scatter(x,v,5,color);
%     [C,h] =contourf(f,50);
%     imagesc(f);
%     colormap(hsv);
    subplot(2,2,1:2);
    colorbar;
    [C,h] =contourf(ax,av,f,100);
    set(h,'Color','none')   
    set(gca,'YDir','normal');
    title('Phase Space, f(x,v,t)');
    text(-2,1,['t','=',num2str(i*dt)]);
    xlabel('x'); ylabel('v');
    
    subplot(2,2,3);
    plot(vv,fv0,'--',vv,fv,'-');
    title('f(v,t)');xlabel('v');ylabel('f');
    ylim([0,1.2*maxfv0]);xlim([-vmax,vmax]);
    subplot(2,2,4);
    plot(vv,fv-fv0,'-',[-vmax,vmax],[0,0],'--');
    title('\delta{}f(v,t)');xlabel('v');ylabel('\delta{}f');
    xlim([-vmax,vmax]);
    ylim([-0.2*maxfv0,0.2*maxfv0]);
    
    F(j)=getframe(gcf,[0,0,800,600]);
    j=j+1;
end
movie(F);
writegif('plot_fxv.gif',F,0.1);
close all;
    
