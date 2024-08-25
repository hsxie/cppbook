% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2016-06-28 11:48
% Compare J0 and N-ring average, with kx=ky=k_perp/sqrt(2)
% Ref: Lee1987, Broemstrup2008 thesis p37
close all; clear; clc;
figure('unit','normalized','position',[0.02,0.2,0.4,0.6],...
    'DefaultAxesFontSize',14);

plt=2;
if(plt==1)
% 1.
M=4;
kp=0:0.1:10;
nk=length(kp);
Jm=0.*kp;
for ik=1:nk
    for jm=1:M
        kx=kp(ik)*cos(jm*2*pi/M)/sqrt(2);
        ky=kp(ik)*sin(jm*2*pi/M)/sqrt(2);
        Jm(ik)=Jm(ik)+exp(1i*(kx+ky))/M;
    end
end
J0=besselj(0,kp);
plot(kp,J0,'-',kp,real(Jm),'r--','Linewidth',2);
legend('J_0',['N=',num2str(M)]); legend('boxoff');
print(gcf,'-dpng',['gk_n_ring_n=',num2str(M),'.png']);
max(abs(imag(Jm))) % if M \neq 2^n, imag(Jm) \neq 0
else
% 2.
kp=0:0.1:30;
nk=length(kp);
calstr=['Jm=0.*kp;for ik=1:nk, for jm=1:M, kx=kp(ik)*cos(jm*2*pi/M)/sqrt(2);',...
    'ky=kp(ik)*sin(jm*2*pi/M)/sqrt(2);Jm(ik)=Jm(ik)+exp(1i*(kx+ky))/M;end,end'];
J0=besselj(0,kp);
M=4;eval(calstr);Jm4=Jm;
M=8;eval(calstr);Jm8=Jm;
M=16;eval(calstr);Jm16=Jm;
M=32;eval(calstr);Jm32=Jm;

plot(kp,J0,'k-',kp,real(Jm4),'g.',kp,real(Jm8),'r--',...
    kp,real(Jm16),'b-.',kp,real(Jm32),'m:','Linewidth',2);
legend('J_0','N=4','N=8','N=16','N=32'); legend('boxoff');
ylim([-1.2,2.4]); xlabel('k_\perp\rho_i'); title('J_0 vs. N-ring average');
print(gcf,'-dpng',['gk_n_ring_cmp.png']);

end