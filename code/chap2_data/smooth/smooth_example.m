close all;clear;clc;
t =0:0.01:10; yt0=sin(pi.*t)+0.5.*cos(2.0*pi.*t);
ytn=yt0+0.5.*(rand (1,length(t))-0.5);
yts1=smooth(ytn,5); % 'moving', 'lowess', 'sgolay', ...
figure; set (gcf,'DefaultAxesFontSize',15);
subplot(211);plot(t,ytn,'LineWidth',2);title('noise data');ylim([-2.5,2.5]);
subplot(212);plot(t,yt0,t,yts1,'r','LineWidth',2);ylim([-2.5,2.5]);
legend('original data','smooth data');legend('boxoff');