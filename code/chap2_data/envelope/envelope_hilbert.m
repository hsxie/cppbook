clear;clc;
t = 0:.05:500;
x0 = sin(t)+1.1.*sin(1.1*t); % signal data
x1=hilbert(x0); % hilbert data, x1=x0+i*xh
plot(t,x0,'g',t,abs(x1),'r');
title('evenlope analysis using Hilbert transform');
legend('original data','evenlope');