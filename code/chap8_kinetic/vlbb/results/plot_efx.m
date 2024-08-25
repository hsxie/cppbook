i=400;
data=load(['efx',num2str(i,'%4.4d')]);
x=data(:,2);
efx=data(:,3);
rhox=data(:,4);
rhox=rhox-mean(rhox);

subplot(211);plot(x,efx);xlim([min(x),max(x)]);
ylim([1.1*min(efx)-0.1*max(efx),1.1*max(efx)-0.1*min(efx)]);
title('E-x');ylabel('E');xlabel('x');

subplot(212);plot(x,rhox);xlim([min(x),max(x)]);
ylim([1.1*min(rhox)-0.1*max(rhox),1.1*max(rhox)-0.1*min(rhox)]);
title('\rho-x');ylabel('\rho');xlabel('x');
