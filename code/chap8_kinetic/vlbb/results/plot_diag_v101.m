% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2011-06-23 05:35
% plot_diag.m, for data analysis of vlasov es 1d code vl1dper.f90
% Time Density Momentum Kinetic Electric Total Entropy1 Entropy2
close all;clear;clc;
h = figure;
diag=load('history.out');
t=diag(:,1);
dens=diag(:,2);
p=diag(:,3);
efmax=diag(:,4);
ek=diag(:,5);
ee=diag(:,6);
et=diag(:,7);
s1=diag(:,8);
s2=diag(:,9);

subplot(221);plotyy(t,dens-mean(dens),t,p-mean(p));xlim([min(t),max(t)]);
title('Density & Momentum');xlabel('t');ylabel('n-n_0 & P-P_0');
legend('n-n_0','P-P_0');

subplot(222);
extrMaxValue = efmax(find(diff(sign(diff(efmax)))==-2)+1);
extrMaxIndex = find(diff(sign(diff(efmax)))==-2)+1;
gamma=-0.153359*pi;efmax0=efmax(1);
% plot(t,log(efmax),t,log(efmax0.*exp(gamma.*t)),'r--');
plot(t,efmax);
% plot(t,efmax,t(extrMaxIndex),efmax(extrMaxIndex),'r*');
% plot(t,log(efmax),t(extrMaxIndex),log(efmax(extrMaxIndex)),'r*');
title('|E|_{max}-Time');xlim([min(t),max(t)]);
xlabel('t');ylabel('|E|_{max}');
Tr=t(extrMaxIndex(2:end))-t(extrMaxIndex(1:end-1));
wr=Tr./(2.0*pi);
wi=(log(efmax(extrMaxIndex(2:end)))-log(efmax(extrMaxIndex(1:end-1))))./Tr;


subplot(223);
% semilogy(t,ek,t,ee,t,et);
% plot(t,log(ek),t,log(ee),t,log(et));
plot(t,ek,t,ee,t,et);xlim([min(t),max(t)]);
title('Energy-Time');xlabel('t');ylabel('E');
legend('E_k','E_e','E_{tot}',1);

subplot(224);
plotyy(t,diag(:,7),t,diag(:,8));xlim([min(t),max(t)]);
title('Entropy-Time');xlabel('t');ylabel('En');
legend('En1','En2');

print(h,'-djpeg','plot_diag');
% close all;

figure;
spectrogram(efmax,128,120,128,1E3);
