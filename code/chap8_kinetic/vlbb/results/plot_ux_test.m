close all;clear;clc;

xl=10.0; dt=0.01;
n=64;nplot=10;nt=1000;
dx=xl/(n-1);

figure('position',[50 50 800 600]);
set(gca,'nextplot','replacechildren');
j=1;
for i=0:nplot:(nt-nplot)

    data=load(['ux',num2str(i,'%4.4d')]);
    x=data(:,2);
    ux=data(:,3);

    plot(x,ux);xlim([min(x),max(x)]);
    ylim([0,1.2]);
    title('u-x');ylabel('u');xlabel('x');
    
    F(j)=getframe(gcf,[0,0,800,600]);
    j=j+1;
end
movie(F);
writegif('plot_ux_test.gif',F,0.1);
close all;