close all;clear;clc;

% xl=10.0; dt=0.01;
n=64;nplot=400;nt=5000;
% dx=xl/(n-1);

figure('position',[50 50 800 600]);
set(gca,'nextplot','replacechildren');
j=1;
for i=0:nplot:(nt-nplot)

    data=load(['efx',num2str(i,'%4.4d')]);
    x=data(:,2);
    efx=data(:,3);
    efx=efx-mean(efx);
    rhox=data(:,4);
    rhox=rhox-mean(rhox);

    subplot(211);plot(x,efx);xlim([min(x),max(x)]);
%     ylim([0,1.2]);
    title('ef-x');ylabel('ef');xlabel('x');
    
    subplot(212);plot(x,rhox);xlim([min(x),max(x)]);
%     ylim([0,1.2]);
    title('\rho-x');ylabel('\rho');xlabel('x');
    
    F(j)=getframe(gcf,[0,0,800,600]);
    j=j+1;
end
movie(F);
writegif('plot_efx_test.gif',F,0.1);
close all;