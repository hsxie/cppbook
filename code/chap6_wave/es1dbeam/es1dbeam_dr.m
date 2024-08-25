% beam-plasma, Hu2006 p465, 1=1/x^2+nb/(x-k*u)^2
close all;clear;clc;
nb=0.1; u0=1; w=[]; kk=-4:0.05:4;
for k=kk
    p=[1, -2*k*u0, (k^2*u0^2-nb-1), 2*k*u0, -k^2*u0^2];
    omg=roots(p); w=[w,omg];
end
wr1=real(w(1,:));wr2=real(w(2,:));wr3=real(w(3,:));wr4=real(w(4,:));
wi1=imag(w(1,:));wi2=imag(w(2,:));wi3=imag(w(3,:));wi4=imag(w(4,:));
figure('DefaultAxesFontSize',15);
plot(kk,wr1,'r.',kk,wi1,'g.',kk,wr2,'r.',kk,wi2,'g.',...
    kk,wr3,'r.',kk,wi3,'g.',kk,wr4,'r.',kk,wi4,'g.');
title(['Beam-plasma dispersion relation, nb=',num2str(nb),...
    ', u0=',num2str(u0)]);
xlabel('k'); ylabel('\omega');
legend('\omega_r','\omega_i',2);
