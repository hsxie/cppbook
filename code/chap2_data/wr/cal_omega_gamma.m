% Hua-sheng XIE, 2015-05-10
close all;clear;clc;
w=3.54+0.23i; nt=1000; dt=0.01; tt=linspace(0,nt*dt,nt+1);
yt=0.1*sin(real(w).*tt).*exp(imag(w).*tt);
figure('DefaultAxesFontSize',15);
subplot(121); plot(tt,yt,'LineWidth',2); hold on;
xlabel('t'); ylabel('yt');
title(['(a) \omega^T=',num2str(real(w)),', \gamma^T=',num2str(imag(w))]);
% Find the corresponding indexes of the extreme max values 
it0=floor(nt*2/10); it1=floor(nt*9/10);
lnyt=log(yt); yy=lnyt(it0:it1);
extrMaxIndex = find(diff(sign(diff(yy)))==-2)+1;
t1=tt(it0+extrMaxIndex(1));t2=tt(it0+extrMaxIndex(end));
y1=yy(extrMaxIndex(1));y2=yy(extrMaxIndex(end));
subplot(122); plot(tt,lnyt,'b',[t1,t2],[y1,y2],'r*--','LineWidth',2);
Nw=length(extrMaxIndex)-1; omega=pi/((t2-t1)/Nw);
gammas=(real(y2)-real(y1))/(t2-t1);
title(['(b) \omega^S=',num2str(omega),', \gamma^S=',num2str(gammas)]);
xlabel('t'); ylabel(['ln(yt), N_w=',num2str(Nw)]);
