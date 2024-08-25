close all;clear;clc;
x=0:0.2*pi:2.0*pi; y1=2.*sin(x)+(rand(1,length(x))-0.5);
y2=2.*sin(x)+(rand(1,length(x))-0.5);
y3=2.*sin(x)+(rand(1,length(x),1)-0.5);
yd=[y1;y2;y3]; y=mean(yd); e=std(yd,1,1);
figure;set(gcf,'DefaultAxesFontSize',15);
errorbar(x,y,e,':bs','LineWidth',2); xlabel('x'); ylabel('y');
hold on; plot(x,y1,'*',x,y2,'*',x,y3,'*');
title('error bar');axis tight;