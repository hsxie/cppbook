close all;clear;clc;
load qpc
nfft=64;
figure; set (gcf,'DefaultAxesFontSize',15);
subplot(221);surfc(zmat);axis tight;grid on;
title('Original data');xlabel('time');ylabel('record #');
ftspec=fft2(zmat, nfft, nfft);
ftspec=fftshift(ftspec);
waxis = [-nfft/2:(nfft/2-1)]'/nfft;
subplot(222);contour(waxis,waxis,ftspec);grid on;
title('FFT');xlabel('record');ylabel('f');

subplot(223);dbspec=bicoher(zmat,nfft);
title('Bicoherence');
subplot(224);dbspec=bicoher(zmat(1:20,:),nfft);
title('Bicoherence, less data');