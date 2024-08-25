Fs = 1e3;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 5e3;                     % Length of signal
t = (0:L-1)*T;                % Time vector
% Sum of a 60*t Hz sinusoid and a 120 Hz sinusoid
x = 0.7*sin(2*pi*60*t.*t) + cos(2*pi*50*t).*exp(2*pi*50.*t*0.001);
y = x + 2*randn(size(t))*1e-1;     % Sinusoids plus noise

h=figure('unit','normalized','position',[0.02,0.1,0.6,0.45],...
    'DefaultAxesFontSize',12);
subplot(121); plot(Fs*t,y); axis tight;
title('Signal Corrupted with Random Noise')
xlabel('Time (ms)')

subplot(122);
[S,F,T,P]=spectrogram(y,128,120,128,1E3);
surf(T,F,10*log10(P),'edgecolor','none'); axis tight; view(0,90);
xlabel('Time (s)'); ylabel('Hz');title('Frequency(t)');