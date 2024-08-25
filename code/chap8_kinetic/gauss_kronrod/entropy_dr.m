% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2016-04-27 18:20
% Solve the local kinetic entropy mode dispersion relation
% Adaptive Simpson quadrature to calculate the 2D integral I_nm

close all; clear; clc;
format long;
runtime=cputime;

% tau=1.0; 
epsn=0.95; % Ln/R
% etai=0;
% wsi=1.0;
% wdi=wsi*epsn;

% fdr=@(w,k) 2.0+(fun_Inm(-w/wdi,k,1,0)*(w-wsi))/wdi+...
%     (fun_Inm(-w/wde,k,1,0)*(w+wsi))/wdi;

kk=0.1:0.1:6.0; x0=0.001+0.003i; 
% kk=2.5:0.1:3.0; x0=0.007+0.05i; 

wk=0.*kk; nk=length(kk);

for jk=1:nk
    
    w=0; k=kk(jk)
    wsi=k; % *v_ti/L_n
    wdi=2*wsi*epsn;
    wse=-wsi;
    wde=-wdi;
    ki=k;
    ke=k/sqrt(1836);
%     ke=k/sqrt(1e6);

    fdr=@(w) 2.0+(fun_gz_gk_Inm(-w/wdi,0,ki^2,1,0)*(w-wsi))/wdi+...
        (fun_gz_gk_Inm(-w/wde,0,ke^2,1,0)*(w-wse))/wde;
      
    options=optimset('Display','off');
    x=fsolve(fdr,x0,options)
    x0=x;
    wk(jk)=x;
end
runtime=cputime-runtime;
kk0=kk; wk0=wk;
%%
close all;
figure('unit','normalized','position',[0.01 0.1 0.6 0.45],...
    'DefaultAxesFontSize',15);
kk=sqrt(2)*kk0;
% ind=find(imag(wk)>1e-6);
% kk=kk(ind); wk=wk(ind);

subplot(121);
plot(kk,real(wk),'bo--','Linewidth',2);
xlabel('k'); ylabel('\omega_r L_n/v_{ti}');
ylim([-0.00,0.006]);
xlim([0,max(kk)]);
title(['entropy mode, L_n/R=',num2str(epsn),', \tau=',...
    num2str(1.0)]);
subplot(122);
plot(kk,imag(wk),'ro--','Linewidth',2);
ylim([-0.0,0.08]); 
xlim([0,max(kk)]);
xlabel('k'); ylabel('\omega_i L_n/v_{ti}');
title(['nk=',num2str(nk),...
    ', runtime=',num2str(runtime),'s']);
% 
print(gcf,'-dpng',['entropy_dr_nk=',...
    num2str(nk),',Ln=',num2str(epsn),'.png']);

