% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2012-01-10 14:01
% On Numerical Calculation of the Plasma Dispersion Function
% Ion acoustic wave, ES1D
clear;clc;
Ti=1; % Ti/Te
mi=1836; % mi/me
zeta=@(x)faddeeva(x)*1i*sqrt(pi);
f=@(x,k)k*k+(1+x*zeta(x))+1/Ti*(1+(x/sqrt(Ti/mi))*zeta(x/sqrt(Ti/mi)));
w=[];

kmin=0.05;dk=0.05;kmax=1.1;
k=kmax:-dk:kmin;
% x0=(1.0-0.0i)*sqrt(2)*kmin;
% x0=1.2-0.5i;
x0=0.04-0.02i;
for kk=k
    options=optimset('Display','off');
    x=fsolve(f,x0,options,kk)*sqrt(2)*kk;
    x0=x/(sqrt(2)*kk);
    w=[w,x];
end
wre=real(w);wie=imag(w);
figure;set(gcf,'DefaultAxesFontSize',15);
plot(k,wre,'r+',k,wie,'.','LineWidth',2);
grid on;
xlabel('k\lambda_D');ylabel('\omega/\omega_p');
xlim([kmin,kmax]); title('ES1D IAW');
legend('\omega_r','\gamma',2); legend('boxoff');