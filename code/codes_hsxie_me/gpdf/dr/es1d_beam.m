% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2012-01-10 14:01
% On Numerical Calculation of the Plasma Dispersion Function, section 10.3
% http://ifts.zju.edu.cn/forum/viewtopic.php?f=18&t=527,bump-on-tail
close all; clear; clc;
nb=0.1; % nb/n0
Tb=1; % Tb/Te
vd=5; % vd/vtb
zeta=@(x)faddeeva(x)*1i*sqrt(pi);
f=@(x,k)k*k+(1-nb)*(1+x*zeta(x))+nb/Tb*(1+(x/sqrt(Tb)-vd)*zeta(x/sqrt(Tb)-vd));
w=[];

kmin=0.05;dk=0.005;kmax=0.4;
% k=[kmin:dk:kmax];
k=kmax:-dk:kmin;
x0=(1.0-0.01i)*sqrt(2)*kmax;
for kk=k
    options=optimset('Display','off');
    x=fsolve(f,x0,options,kk)*sqrt(2)*kk;
    x0=x/(sqrt(2)*kk);
%     x0=x;
    w=[w,x];
end
wre=real(w);wie=imag(w);
figure;
set(gcf,'DefaultAxesFontSize',15);
% subplot(211);
plot(k,wie,'.','LineWidth',2);
grid on;
title(['Bump-on-tail, n_b=',num2str(nb),', T_b/T_e=',...
    num2str(Tb),', v_d=',num2str(vd)]);
xlabel('k\lambda_D');ylabel('\gamma/\omega_p'); % wp^2=wpe^2+wpb^2
xlim([kmin,kmax]);
ylim([-0.08,0.28]);
% subplot(212);
% plot(k,wre,'.','LineWidth',2);


print('-dpng',['beam_nb',num2str(nb),'.png']);