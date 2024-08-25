% Hua-sheng XIE, huashengxie@gmail.com, FSC-ZJU, 2016-04-29 17:30
% matrix to solve local kinetic ITG dispersion relation
% lambda*[g_ij(v_par,v_perp),phi]=M_ij*[g_ij(v_par,v_perp),phi]
% Benchmark with Kim1994PoP Fig.5
% 16-05-10 09:21 Fix bug of F0
% 16-10-17 13:02 rewrite for entropy mode
% 17-01-13 10:51 tau\neq1.0, bug on phi=sum(gi+ge) should be sum(gi+ge/tau),
% 11:36 fixed, test tau=2.0 agree with IVP & DR
close all;clear;%clc;

data =[0.1,     0.1277+2.8258i;
       0.2,     0.2487+2.8211i;
       0.5,     0.6369+2.6664i;
       1.0,     1.3108+2.2697i;
       1.5,     2.1182+1.9633i;
       2.0,     3.1047+1.9143i];
id=4;
ky=data(id,1); % k_perp
wr=real(data(id,2)); wi=imag(data(id,2));

runtime=cputime;

tau=1.0;
epsn=0.2;
kapn=1/epsn; % R/L_n
kapt=0.1*kapn; % R/L_T
eta=kapt/kapn;

kz=0.0; % k_z*R
mi=1836.0; % mi/me

wdi=ky; % *v_ti/R
wde=-wdi*tau;
ki=ky;
ke=-ky*sqrt(tau/mi);
kzi=kz;
kze=kz*sqrt(tau*mi);
Gamma0i=besseli(0,ki^2,1);
Gamma0e=besseli(0,ke^2,1);
G0coef=1.0/((1.0-Gamma0i)+(1.0-Gamma0e)/tau)/sqrt(2.0*pi);

nvx=1*64/1; nvy=1*64/1;
vxmax=5.0; vxmin=-vxmax; vymax=5.0; vymin=0.0;
% vxmax=6.0; vxmin=-vxmax; vymax=6.0; vymin=0.0;
% vxmax=10.0; vxmin=-vxmax; vymax=10.0; vymin=0.0;

dvx=(vxmax-vxmin)/nvx; dvy=(vymax-vymin)/nvy;

Vx=vxmin:dvx:vxmax; Vy=vymin:dvy:vymax; 
[vx,vy]=meshgrid(Vx,Vy); % vx and vy grids
wDiv=wdi*(vx.^2+vy.^2/2);
wTiv=wdi*(kapn+(0.5*(vx.^2+vy.^2)-1.5)*kapt);
wDev=wde*(vx.^2+vy.^2/2);
wTev=wde*(kapn+(0.5*(vx.^2+vy.^2)-1.5)*kapt);

F0=exp(-0.5*(vx.^2+vy.^2)); % initial distribution

[Nx,Ny]=size(vx);
Nda=Nx*Ny;
Nd=2*Nx*Ny+1; II=[]; JJ=[]; SS=[];
cout=0;
for jx=1:Nx
    for jy=1:Ny
        
        % for ion
        ind=(jx-1)*Ny+jy; 
        
        cout=cout+1;
        II(cout)=ind; JJ(cout)=ind;
        SS(cout)=(kzi*vx(jx,jy)+wDiv(jx,jy));
        
        cout=cout+1;
        II(cout)=ind; JJ(cout)=Nd;
        SS(cout)=((kzi*vx(jx,jy)+wDiv(jx,jy)-wTiv(jx,jy))*(...
            besselj(0,ki*vy(jx,jy)))^2*F0(jx,jy));
        
        cout=cout+1;
        II(cout)=Nd; JJ(cout)=ind;
        SS(cout)=G0coef*(kzi*vx(jx,jy)+wDiv(jx,jy))*vy(jx,jy)*dvx*dvy;
        
        % for electron
        ind=ind+Nda; 
        
        cout=cout+1;
        II(cout)=ind; JJ(cout)=ind;
        SS(cout)=(kze*vx(jx,jy)+wDev(jx,jy));
        
        cout=cout+1;
        II(cout)=ind; JJ(cout)=Nd;
        SS(cout)=((kze*vx(jx,jy)+wDev(jx,jy)-wTev(jx,jy))*(...
            besselj(0,ke*vy(jx,jy)))^2*F0(jx,jy));
        
        cout=cout+1;
        II(cout)=Nd; JJ(cout)=ind;
        SS(cout)=G0coef/tau*(kze*vx(jx,jy)+wDev(jx,jy))*vy(jx,jy)*dvx*dvy; % 17-01-13 11:29 fixed bug on tau
        
    end
end
cout=cout+1;
II(cout)=Nd; JJ(cout)=Nd;
SS(cout)=G0coef*sum(sum(((kzi*vx+wDiv-wTiv).*(...
        besselj(0,ki.*vy)).^2+(kze*vx+wDev-wTev).*(...
        besselj(0,ke.*vy)).^2/tau).*F0.*vy))*dvx*dvy; % 17-01-13 10:52 fixed bug on tau

M=sparse(II,JJ,SS,Nd,Nd);
% wg=1.2 + 2.3i; % k=0.9
wg=wr+1i*wi;
if (nvx*nvy<1000)
    d=eig(full(M));
else
    d=eigs(M,1,wg);
end
% d=eigs(M,1,'li');
runtime=cputime-runtime;
%%
[di,ind]=sort(imag(d),'descend');
d=d(ind);
% d(1:5) % find(imag(d)>0)

[V,D]=eigs(M,1,d(1));
w=D(1,1);
X=V(:,1);
%%
gi=reshape(V(1:Nda),Ny,Nx);
ge=reshape(V((Nda+1):(Nd-1)),Ny,Nx);
gi=gi.'; ge=ge.';

close all;
h = figure('Unit','Normalized','position',...
    [0.01 0.13 0.6 0.65],'DefaultAxesFontSize',15);

subplot(221);
plot(real(d),imag(d),'x',real(w),imag(w),'rs','Linewidth',2);
xlabel('\omega_r'); ylabel('\omega_i');
title(['(a) k_{||}=',num2str(kz),', k_\perp=',num2str(ky),...
    ', \epsilon_n=',num2str(epsn),', \eta=',num2str(eta),...
    ', \tau=',num2str(tau)]);

subplot(222);
% pcolor(vx,vy,real(ge)); shading interp; %log(abs(real(g)))
contourf(vx,vy,real(ge),30,'LineStyle','none'); colorbar;
xlabel('v_{||}'); ylabel('v_\perp'); 
title(['(b) Re g_e, \omega=',num2str(w,4)]);
xlabel(['\omega^T=',num2str(wr+1i*wi,4)]);

subplot(223);
% pcolor(vx,vy,real(g)); shading interp; %log(abs(real(g)))
contourf(vx,vy,real(gi),30,'LineStyle','none'); colorbar;
xlabel(['v_{||}, runtime=',num2str(runtime),'s']); ylabel('v_\perp'); 
title(['(c) Re g_i, nvx=',num2str(nvx),', nvy=',num2str(nvy)]);
subplot(224);
% pcolor(vx,vy,imag(g)); shading interp;
contourf(vx,vy,imag(gi),30,'LineStyle','none'); colorbar;
xlabel('v_{||}'); ylabel('v_\perp'); 
title(['(d) Im g_i, vxmax=',num2str(vxmax),', vymax=',num2str(vymax)]);

str=['entropy_matrix_kz=',num2str(kz),',ky=',num2str(ky),...
    ',epsn=',num2str(epsn),',eta=',num2str(eta),...
    ',tau=',num2str(tau),',vxmax=',num2str(vxmax),...
    ',vymax=',num2str(vymax),',nvx=',num2str(nvx),...
    ',nvy=',num2str(nvy)];

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6.0]);
print(gcf,'-dpng',[str,'.png'],'-r100');
% print(gcf,'-dpdf',[str,'.pdf'],'-r100');
