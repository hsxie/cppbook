close all;clear;clc;
h = figure;
diag=load('history.out');
t=diag(:,1);
efmax=diag(:,4);

subplot(211);
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

subplot(212);
spectrogram(efmax,256*2);

print(h,'-djpeg','plot_diag');
% close all;

