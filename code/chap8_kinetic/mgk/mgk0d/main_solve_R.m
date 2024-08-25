% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2016-10-10 12:06
% Solve the local kinetic entropy mode dispersion relation
% Adaptive Simpson quadrature to calculate the 2D integral
% Rewrite from the 2016-04-27 18:20 version
% 16-10-11 09:52 change the output normalization, with R/vti not L_n/vti
% 16-10-11 14:20 this version fixed the possible bug at small k, by add
% the tolerace to 1e-6 to 1e-10 in fsolve()
% 16-10-14 16:03 compare z-pinch and dipole
% 16-10-16 16:14 with kpara

close all; clear; clc;
runtime=cputime;

tol=1e-6; xmax=1.0e1; ymax=1.0e1; %  control the 2d integral accuracy 

tau=1.0;
epsn=0.2;
kapn=1/epsn; % R/L_n
kapt=0.1*kapn; % R/L_T
eta=kapt/kapn;

kz=0.0;
mi=1836.0; % mi/me

% kk=4.0:-0.2:0.2; x0=0.01+0.03i;
kk=0.1:0.1:2.0; x0=0.01+0.03i;
wk=0.*kk; nk=length(kk);

for jk=1:nk
    
    w=0; k=kk(jk)
    wdi=k; % *v_ti/R
    wde=-wdi*tau;
    ki=k;
    ke=-k*sqrt(tau/mi);
    kzi=kz;
    kze=kz*sqrt(tau*mi);
    
    fdr=@(w) fun_entropy_dr(w,ki,wdi,kapn,kapt,kzi,tol,xmax,ymax)+...
        (1/tau)*fun_entropy_dr(w,ke,wde,kapn,kapt,kze,tol,xmax,ymax);
      
    options=optimset('Display','off');
%     options=optimset('Display','off','TolX',1e-10,'TolFun',1e-10); % 16-10-11 14:20
    x=fsolve(fdr,x0,options)
    x0=x;
    wk(jk)=x;
end
runtime=cputime-runtime;
%%
dat0=[kk.',wk.']; wk0=wk;
dat=dat0(1:2:end,:);
%%
ind=find(imag(wk)<1e-4 | abs(real(wk))>10e0); wk(ind)=NaN+1i*NaN;
close all;
figure('unit','normalized','position',[0.01 0.1 0.6 0.4],...
    'DefaultAxesFontSize',15);

str=['epsn=',num2str(epsn),',eta=',num2str(eta),...
    ',tau=',num2str(tau),',k_z=',num2str(kz),',nk=',num2str(nk)];
str2=[',xmax=',num2str(xmax),',ymax=',num2str(ymax),',tol=',num2str(tol)];

subplot(121);
plot(kk,real(wk),'bd--','Linewidth',2);
xlabel(['k\rho_i', str2]); ylabel('\omega_r R/v_{ti}');
xlim([0,max(kk)]);
title(['entropy mode',',runtime=',num2str(runtime,4),'s']);
% legend('k_z=0',['k_z=',num2str(kz)],1);legend('boxoff');
% ylim([0,0.1]);
subplot(122);
plot(kk,imag(wk),'bd--','Linewidth',2);
xlim([0,max(kk)]);
% ylim([0,0.3]);
xlabel('k\rho_i'); ylabel('\omega_i R/v_{ti}');
title(str);
dat=[kk.',wk.']
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 3.0]);
print(gcf,'-dpng',['entropy_dr_',str,str2,'.png'],'-r100');
print(gcf,'-dpdf',['entropy_dr_',str,str2,'.pdf'],'-r100');


