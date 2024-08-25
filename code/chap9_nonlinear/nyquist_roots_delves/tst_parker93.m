% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-10-08 17:19
% Using Delves-Lyness 1967 approach to solve Parker93 D.R.
% Define the f(z)=0 function at funf.m

close all; clear; clc;

global mi Ti vti me tau Te vte alphai alphae theta k kappan kappat;

mi=1.0; Ti=1.0; vti=sqrt(Ti/mi);
me=1/1837; tau=1.0; Te=tau*Ti; vte=sqrt(Te/me);
alphai=1; alphae=-1/me; theta=0.01;
k=0.5*2; kappan=0.0; kappat=2.0;

runtime=cputime;

M=1;
% za=-0.1-0.03i; zb=0.7+0.2i; % paker93
% za=-0.1-0.03i; zb=0.2+0.02i; % paker93
za=-6.1-4.0i; zb=8.2-0.0i; % landau
% za=-10.1-5.0i; zb=13.2-0.0i; % landau
% za=-0.5+0.01i; zb=2.0+1.3i; % paker93
% za=-0.5+0.01i; zb=1.1/2+0.5i; % paker93
% za=-1.51-0.01i; zb=1.31+1.1i; % paker93

xa=real(za); ya=imag(za); xb=real(zb); yb=imag(zb);

[N,domain]=fun_divide_domain(za,zb,M);
runtime=cputime-runtime;

%%
close all;
h=figure('unit','normalized','defaultaxesfontsize',15,...
    'position',[0.02 0.1 0.5 0.55]);

ndom=length(domain);
ww=[];
for jd=1:ndom
    zda=domain(jd).za; zdb=domain(jd).zb;
    tol=1e-4;
    [Nd,rtd]=fun_rt(zda,zdb,tol);
    ww=[ww;rtd];
    plot(real(rtd),imag(rtd),'rx','Markersize',5,'linewidth',2); hold on;
    xda=real(zda); xdb=real(zdb);
    yda=imag(zda); ydb=imag(zdb);
    xpd=[xda,xdb,xdb,xda,xda];
    ypd=[yda,yda,ydb,ydb,yda];
    plot(xpd,ypd,'c:','linewidth',2); hold on;
    text(real(rtd(1))+0.1,imag(rtd(1)),num2str(rtd(1),4));
end
[wi,jw]=sort(imag(ww),'descend');
ww=ww(jw);
plot([xa,xb,xb,xa,xa],[ya,ya,yb,yb,ya],'--','linewidth',2); hold on;
xlabel('\omega_r'); ylabel('\omega_i');

str=['M=',num2str(M),',N=',num2str(N),...
    ',za=',num2str(za),',zb=',num2str(zb)];
title([str,', runtime=',num2str(runtime),'s']);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 4.5]);
print(gcf,'-dpng',['divide_domain_Landau_roots_',str,'.png'],'-r100');
% print(gcf,'-dpng',['divide_domain_Parker93_roots_',str,'.png'],'-r100');
