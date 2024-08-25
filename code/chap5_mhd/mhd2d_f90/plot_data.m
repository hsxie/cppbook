% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2012-06-17 11:54
% plot_data.m, for data analysis of mhd 2d code mhd2d.f90
close all;clear;clc;

% modify here base on parameters of mhd2d.f90
dx=1.0; dz=2.0; dt=0.05;
ni=2^5+1; nj=2^5; nt=5000; nplot=floor((nt+1)/50);

xl=dx*(ni-1)/2; zl=dz*(nj-1);
xx=-xl:dx:xl;
[ax,az]=meshgrid(-xl:dx:xl,0:dz:zl);
data0=load('bxm.dat');
itt=data0(:,1);
tt=data0(:,2);
bxmax=data0(:,3);

figure('position',[50 50 800 600]);
set(gca,'nextplot','replacechildren');
set(gcf,'DefaultAxesFontSize',15);
j=1;
for i=0:nplot:(nt-nplot)
    
    data1=load(['rho',num2str(i,'%6.6d')]);
    data2=load(['ux',num2str(i,'%6.6d')]);
    data3=load(['Bx',num2str(i,'%6.6d')]);
    x=data1(:,3);
    z=data1(:,4);
    rho=reshape(data1(:,5),nj,[]);
    ux=reshape(data2(:,5),nj,[]);
    Bx=reshape(data3(:,5),nj,[]);
    
    subplot(221); pcolor(ax,az,rho); shading('interp'); % pcolor
    title(['rho(x,z,t)',', t=',num2str(i*dt)]);
    xlabel('x'); ylabel('z');
    
    subplot(222); [C,h] =contourf(ax,az,ux,100);
    set(h,'Color','none'); set(gca,'YDir','normal');
    title(['ux(x,z,t)',', t=',num2str(i*dt)]);
    xlabel('x'); ylabel('z');
    
    subplot(223); [C,h] =contourf(ax,az,Bx,100);  % contourf
    set(h,'Color','none'); set(gca,'YDir','normal');
    title(['Bx(x,z,t)',', t=',num2str(i*dt)]);
    xlabel('x'); ylabel('z');
    
    subplot(224);
    plot(tt,log(bxmax),'g:',tt(1:i),log(bxmax(1:i)),'r','LineWidth',2);
    title(['max(Bx(t))',', t=',num2str(i*dt)]);
    xlabel('time'); ylabel('log(max(Bx))'); xlim([0,max(tt)]);

    F(j)=getframe(gcf,[0,0,800,600]);
    j=j+1;
end
movie(F);
writegif('plot_data.gif',F,0.1);
close all;
    
