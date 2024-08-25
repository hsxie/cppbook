% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2016-05-09 17:17
% Solve the local kinetic ITG mode dispersion relation
% Adaptive Simpson quadrature to calculate the 2D integral
% 16-09-30 08:30 update

close all; clear; clc;
runtime=cputime;

kapn=1.0;
kapt=2.5;
epsn=0.2;

% kk=0.05:0.05:1.2; x0=0.2+0.3i;
kk=0.1:0.1:1.2; x0=0.02+0.06i;
% kk=0.5:-0.05:0.05; x0=0.2+0.3i;
% kk=0.3:0.05:2.0; x0=-0.2+0.6i;
% kk=[0.22,0.25:0.05:1.0,1.1:0.1:1.8]; x0=-0.02+0.006i;
% kk=2.0:0.05:3.0; x0=0.64+0.039i;
wk=0.*kk; nk=length(kk);
for jk=1:nk
    
    k=kk(jk)
    ky=k;
    ws=ky; % *v_ti/L_n
    wd=ky*epsn;
%     wd=0;
        
    fdr=@(w) 2.0-fun_as_itgdr(w,ky,kapn,kapt,ws,wd);
      
    options=optimset('Display','off');
    x=fsolve(fdr,x0,options)
    x0=x;
    wk(jk)=x;
end
runtime=cputime-runtime;
kk0=kk;
%%
% kk0=kk;
% kk=kk0/sqrt(2);

ind=find(wk~=0);
kk=kk(ind);
wk=wk(ind);
dat=[kk',(real(wk))',(imag(wk))'];

% close all;
figure('unit','normalized','position',[0.01 0.1 0.6 0.45],...
    'DefaultAxesFontSize',15);

subplot(121);
plot(kk,real(wk),'bo--','Linewidth',2);
xlabel('k'); ylabel('\omega_r L_n/v_{ti}');
% ylim([-0.1,0.8]); 
xlim([0,max(kk)]);
title(['itg, L_n/R=',num2str(epsn),', \tau=',...
    num2str(1.0),', \kappa_n=',num2str(kapn),', \kappa_t=',num2str(kapt)]);
subplot(122);
plot(kk,imag(wk),'ro--','Linewidth',2);
% ylim([-0.01,0.2]);
xlim([0,max(kk)]);
xlabel('k'); ylabel('\omega_i R/v_{ti}');
title(['nk=',num2str(nk),...
    ', runtime=',num2str(runtime),'s']);

print(gcf,'-dpng',['itg_dr_epsn=',num2str(epsn),',kapn=',num2str(kapn),...
    ',kapt=',num2str(kapt),',nk=',num2str(nk),'_2.png']);
