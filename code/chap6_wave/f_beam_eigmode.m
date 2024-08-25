% Hua-sheng XIE, 17-01-05 16:13, f_beam_eigmode.m
% To show the existence of three types of spectral (continuum, discrete and residual modes)
% in Vlasov-Ampere system with bump-on-tail distribution function 
close all; clear; clc;

k=0.2;vmax=8.0; N=128*2/1;
v = linspace(-vmax,vmax,N+1); % row of the vector v
dv = v(2) - v(1);

%%
Nv=N+1;
Mva=zeros(Nv+1);

nb=0.2; vtm=1.0; vtb=1*vtm; vb=5*vtm;% Maxwellian or bump-on-tail
for j=1:Nv
    f0(j)=(1-nb)/vtm/sqrt(2*pi)*exp(-0.5*(v(j)/vtm)^2)+...
        nb/vtb/sqrt(2*pi)*exp(-0.5*((v(j)-vb)/vtb)^2);
    dvf0(j)=-v(j)*(1-nb)/vtm^3/sqrt(2*pi)*exp(-0.5*(v(j)/vtm)^2)-...
        (v(j)-vb)*nb/vtb^3/sqrt(2*pi)*exp(-0.5*((v(j)-vb)/vtb)^2);
end

for j=1:Nv
    Mva(j,j)=k*v(j);
    Mva(j,Nv+1)=1i*dvf0(j);
    Mva(Nv+1,j)=1i*v(j)*dv;
end

w=eig(Mva);
[neww,inds] = sort(real(w));

h=figure('unit','normalized','Position',[0.01 0.07 0.6 0.75]);
set(gcf,'DefaultAxesFontSize',15);

subplot(235);
plot(v,f0,'LineWidth',2);
xlabel('v_j'); ylabel('f_0');

subplot(236);
text(0.2,0.5,['v_{max}=',num2str(vmax),10,'N_v=',num2str(N),10,...
    'n_b=',num2str(nb),10,'v_b=',num2str(vb),10,...
    'v_{tb}=',num2str(vtb),10,'k=',num2str(k)],'Fontsize',12);
axis off;

subplot(231); plot(real(w),imag(w),'x','LineWidth',2); 
xlabel('\omega_r'); ylabel('\omega_i'); hold on;
title('(a) V-A eigenmode solutions');

subplot(232);
[VV0,MM0]=eigs(sparse(Mva),2,'li');
plot(v,real(VV0(1:(end-1),1)),v,imag(VV0(1:(end-1),1)),'--','LineWidth',2);
xlabel('v_j, \gamma>0 discrete mode'); ylabel('\delta f(v_j)');
title(['(c) \omega=',num2str(MM0(1,1),3)]);

[VV,MM]=eigs(sparse(Mva),6,0.01);
subplot(233);
plot(v,real(VV(1:(end-1),2)),v,imag(VV(1:(end-1),2)),'--','LineWidth',2);
xlabel('v_j, residual mode'); ylabel('\delta f(v_j)');
title('(c) \omega=0');
subplot(234);
plot(v,real(VV(1:(end-1),1)),v,imag(VV(1:(end-1),1)),'--','LineWidth',2);
xlabel('v_j, continuum mode'); ylabel('\delta f(v_j)');
hl2=legend('Re[\delta f]','Im[\delta f]',1); 
set(hl2,'Fontsize',10,'Box', 'off');
title('(d) \omega=0 CvK');

%%
% plot contour of the Landau solution
wrmax=2; wimax=0.5; wrmin=-wrmax; wimin=-2*wimax; 
vc = 0:0.001:0.01;
zeta=@(x) faddeeva(x)*1i*sqrt(pi);
zetay=@(x)(1+x.*zeta(x));
[Rew,Imw] = meshgrid(wrmin:.002:wrmax, wimin:.002:wimax);
f1 = k^2 + (1-nb)*zetay((Rew+1i*Imw)/(sqrt(2)*k))+ nb/vtb^2*zetay((Rew+1i*Imw-k*sqrt(1)*vb)/(sqrt(2)*k*vtb));
f = sqrt(real(f1).^2 + imag(f1).^2);
subplot(231); hold on;
contour(Rew, Imw, f, vc);
hl=legend('eigenmodes','Landau solutions',4);
set(hl,'Fontsize',6,'Box','off');
ylim([wimin,wimax]); xlim([-3,4]);

print('-dpng','eigenmode.png');
