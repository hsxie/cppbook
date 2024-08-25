close all;clear;clc;
load('data.mat');
figure('unit','normalized','position',[0.02,0.1,0.6,0.6],...
    'DefaultAxesFontSize',12);
subplot(221); plot(t,a); axis tight;
xlabel('time (ms)'); ylabel('A'); title('(a) Signal');
subplot(222); specgram(a,1024,1e3,hanning(512),120);
xlabel('time (ms)'); ylabel('f (kHz)'); title('(b) Frequency');

subplot(223); plot(t,a); axis tight;
xlabel('time (ms)'); ylabel('A'); title('(c) Signal'); xlim([500,700]); ylim([-0.5,0.5]);
subplot(224); specgram(a,1024,1e3,hanning(512),120);
 xlim([500,700]); ylim([0,40]);
xlabel('time (ms)'); ylabel('f (kHz)'); title('(d) Frequency');
