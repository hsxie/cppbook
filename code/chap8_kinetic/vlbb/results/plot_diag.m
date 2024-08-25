% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2011-06-23 05:35
% plot_diag.m, for data analysis of vlasov es 1d code vl1dper.f90
% Time Density Momentum Kinetic Electric Total Entropy1 Entropy2
close all;clear;clc;
h = figure;
diag=load('history.out');
subplot(221);
t=diag(:,1);
plot(t,diag(:,2)-mean(diag(:,2)));
title('Density-Time');xlabel('t');ylabel('n-n_0');
subplot(222);
plot(t,diag(:,3)-mean(diag(:,3)));
title('Momentum-Time');xlabel('t');ylabel('P-P_0');
subplot(223);
gamma=-8.0*0.153359;ee0=diag(1,5);
% semilogy(t,diag(:,4),t,diag(:,5),t,diag(:,6),t,ee0.*exp(gamma.*t),'r--');
% plot(t,log(diag(:,4)),t,log(diag(:,5)),t,log(diag(:,6)),t,log(ee0.*exp(gamma.*t)),'r--');
plot(t,log(diag(:,4)),t,log(diag(:,5)),t,log(diag(:,6)));
% plot(t,(diag(:,4)),t,(diag(:,5)),t,(diag(:,6)));
% plot(t,diag(:,5));
title('Energy-Time');xlabel('t');ylabel('E');
legend('E_k','E_e','E_{tot}',4);
subplot(224);
plot(t,diag(:,7),t,diag(:,8));
title('Entropy-Time');xlabel('t');ylabel('En');
legend('En1','En2','Location','East');
Ee=diag(:,5);
wi=(Ee(2:end)-Ee(1:end-1))./Ee(1:end-1);
wimean=mean(wi(end/8:end*3/8));
print(h,'-djpeg','plot_diag');
% close all;