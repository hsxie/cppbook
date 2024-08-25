% Hua-sheng XIE, huashengxie@gmail.com, FSC-ZJU, 2016-04-29 16:37
% IVP to solve local kinetic entropy mode dispersion relation, 
% d_t g(v_par,v_perp)
% Benchmark with Ricci2006PoP Fig.4
% 16-05-10 23:34
% 16-10-10 15:37 Redo the normalization, omega -> vti/R, k -> k*rhoi
% 17:15 also fixed a bug in RK4
% 18:40 Test OK for Ricci06 Fig.4, k=0.9
% 16-10-11 16:22 fixed another bug in RK4, gi,ge should be gitmp,getmp
% in 2-4 steps
% 16-10-18 17:06 add kpara and collision

close all;clear;clc;

irun=0;
nt=500;
runtimea=cputime;

close all;
if(irun==0)
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

    nvx=1*32/1; nvy=1*64/1; dt=0.02/1; % parameters
    vxmax=5.0; vxmin=-vxmax; vymax=5.0; vymin=0.0;
%     vxmax=3.0; vxmin=-vxmax; vymax=3.0; vymin=0.0;
    % vymax=5.0; vymin=-vymax;

    dvx=(vxmax-vxmin)/nvx;
    dvy=(vymax-vymin)/nvy;
    % fperp=0.1;

    Vx=vxmin:dvx:vxmax; Vy=vymin:dvy:vymax; 
    [vx,vy]=meshgrid(Vx,Vy); % vx and vy grids
    wDiv=wdi*(vx.^2+vy.^2/2)+kzi*vx; % add kpar 16-10-18 17:11
    wTiv=wdi*(kapn+(0.5*(vx.^2+vy.^2)-1.5)*kapt);
    wDev=wde*(vx.^2+vy.^2/2)+kze*vx; %
    wTev=wde*(kapn+(0.5*(vx.^2+vy.^2)-1.5)*kapt);

    F0=exp(-0.5*(vx.^2+vy.^2)); % initial distribution
    % F0=2/sqrt(pi)*vy.*exp(-0.5*(vx.^2+vy.^2).^2); % initial distribution
    gi=0.001.*F0;
    ge=0.001.*F0;

    phi=0.1;
%     phit=zeros(nt,1);
    rk=1;
    
    nta=1; ntb=nt;
else
    nta=it+1; ntb=nta+nt-1;
end

figure('Unit','Normalized','position',...
    [0.02 0.1 0.6 0.4],'DefaultAxesFontSize',15);
for it=nta:ntb
    if(rk==0)
        % boundary in V space g must be null
    %     g(:,1)=0; g(:,nvy+1)=0; g(1,:)=0; g(nvx+1,:)=0;

        phi=G0coef*sum(sum((gi+ge/tau).*vy))*dvx*dvy;
        gi=gi-1i*(wDiv.*gi+(wDiv-wTiv)*phi.*(...
            besselj(0,ki*vy)).^2.*F0)*dt;
        ge=ge-1i*(wDev.*ge+(wDev-wTev)*phi.*(...
            besselj(0,ke*vy)).^2.*F0)*dt;
    else
        % RK-4, 1st step
        dgi1=-1i*(wDiv.*gi+(wDiv-wTiv)*phi.*(...
            besselj(0,ki*vy)).^2.*F0);
        dge1=-1i*(wDev.*ge+(wDev-wTev)*phi.*(...
            besselj(0,ke*vy)).^2.*F0);
        % RK-4, 2nd step
        gitmp=gi+0.5*dt*dgi1;
        getmp=ge+0.5*dt*dge1;
        phitmp=G0coef*sum(sum((gitmp+getmp/tau).*vy))*dvx*dvy;
     %  dgi2=-1i*(wDiv.*gi+(wDiv-wTiv)*phitmp.*(... % wrong, 16-10-11 16:27
     %       besselj(0,ki*vy)).^2.*F0);
        dgi2=-1i*(wDiv.*gitmp+(wDiv-wTiv)*phitmp.*(...
            besselj(0,ki*vy)).^2.*F0);
        dge2=-1i*(wDev.*getmp+(wDev-wTev)*phitmp.*(...
            besselj(0,ke*vy)).^2.*F0);
        % RK-4, 3rd step
        gitmp=gi+0.5*dt*dgi2;
        getmp=ge+0.5*dt*dge2;
        phitmp=G0coef*sum(sum((gitmp+getmp/tau).*vy))*dvx*dvy;
        dgi3=-1i*(wDiv.*gitmp+(wDiv-wTiv)*phitmp.*(...
            besselj(0,ki*vy)).^2.*F0);
        dge3=-1i*(wDev.*getmp+(wDev-wTev)*phitmp.*(...
            besselj(0,ke*vy)).^2.*F0);
        % RK-4, 4th step
%         gitmp=gi+0.5*dt*dgi3; % wrong, 16-10-10 17:15
%         getmp=ge+0.5*dt*dge3;
        gitmp=gi+dt*dgi3;
        getmp=ge+dt*dge3;
        phitmp=G0coef*sum(sum((gitmp+getmp/tau).*vy))*dvx*dvy;
        dgi4=-1i*(wDiv.*gitmp+(wDiv-wTiv)*phitmp.*(...
            besselj(0,ki*vy)).^2.*F0);
        dge4=-1i*(wDev.*getmp+(wDev-wTev)*phitmp.*(...
            besselj(0,ke*vy)).^2.*F0);
        % RK-4, push
        gi=gi+dt./6.0.*(dgi1+2.0*dgi2+2.0*dgi3+dgi4);
        ge=ge+dt./6.0.*(dge1+2.0*dge2+2.0*dge3+dge4);

        phi=G0coef*sum(sum((gi+ge/tau).*vy))*dvx*dvy;
    end
    phit(it)=phi;
    
    if(mod(it,10)==9)
        subplot(121);
        plot((1:it)*dt,real(phit),(1:it)*dt,imag(phit),'r--','LineWidth',2); box on;
        xlabel('t'); legend('Re\phi','Im\phi',2); legend('boxoff');
%         xlim([0,it*dt]);
        title(['\phi, it=',num2str(it),'/',num2str(ntb)]);
        subplot(122);
        plot((1:it)*dt,log(real((phit))),'LineWidth',2); box on;
        xlabel('t'); ylabel('log(phit)');

        drawnow;
    end
end
runtime=runtime+cputime-runtimea;

%%
close all;
h = figure('Unit','Normalized','position',...
    [0.02 0.13 0.6 0.7],'DefaultAxesFontSize',15);

t=(1:ntb)*dt;

subplot(221);
plot(t,real(phit),t,imag(phit),'r--','LineWidth',2); box on;
xlabel(['t',', runtime=',num2str(runtime,4),'s']); 
legend('Re\phi','Im\phi',2); legend('boxoff');
title(['(a) k_\perp=',num2str(ky),', k_z=',num2str(kz),...
    ', \epsilon_n=',num2str(epsn),', \eta_i=',num2str(eta),...
    ', \tau=',num2str(tau)]);

% Find the corresponding indexes of the extreme max values
% lndE=log(sqrt(real((phit(1:nt))))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
lndE=log(real((phit(1:ntb))));
lndEi=log(imag((phit(1:ntb))));
lndEa=log(abs((phit(1:ntb))));
it0=floor(ntb*5.8/20); it1=floor(ntb*19.9/20);
yy=lndE(it0:it1);
subplot(222);
plot(t,lndE,t,lndEi,'b--',t,lndEa,'k:','LineWidth',2); box on;
hold on; axis tight;
withpeak=1;
if(withpeak==1)
    extrMaxIndex = find(diff(sign(diff(yy)))==-2)+1;
    t1=t(it0+extrMaxIndex(1));t2=t(it0+extrMaxIndex(end));
    y1=yy(extrMaxIndex(1));y2=yy(extrMaxIndex(end));
    nw=length(extrMaxIndex)-1;
    omega=pi/((t2-t1)/nw);
else
    t1=t(it0);t2=t(it1);
    y1=yy(1);y2=yy(end);
    nw=0;
    omega=0;
end
plot([t1,t2],[y1,y2],'r*--','LineWidth',2);
gammas=(real(y2)-real(y1))/(t2-t1);
title(['(b) \omega^S=',num2str(omega),', \gamma^S=',num2str(gammas),...
    ', nw=',num2str(nw)]);
xlabel(['\omega^T=',num2str(wr+1i*wi)]);

subplot(223);
% pcolor(vx,vy,real(g)); shading interp; %log(abs(real(g)))
contourf(vx,vy,real(gi),30,'LineStyle','none');
xlabel('v_{||}'); ylabel('v_\perp'); 
title(['(c) Re g_i, nvx=',num2str(nvx),', nvy=',num2str(nvy)]);
subplot(224);
% pcolor(vx,vy,imag(g)); shading interp;
contourf(vx,vy,imag(gi),30,'LineStyle','none');
xlabel('v_{||}'); ylabel('v_\perp'); 
title(['(d) Im g_i, vxmax=',num2str(vxmax),', vymax=',num2str(vymax)]);

str=['mgk0d_ivp_ky=',num2str(ky),...
    ',epsn=',num2str(epsn),',eta=',num2str(eta),...
    ',tau=',num2str(tau),',kz=',num2str(kz),',vxmax=',num2str(vxmax),...
    ',vymax=',num2str(vymax),',nvx=',num2str(nvx),...
    ',nvy=',num2str(nvy),',dt=',num2str(dt)];
print(gcf,'-dpng',[str,'.png']);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6.0]);
print(gcf,'-dpng',[str,'.png'],'-r100');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7.5 5.0]);
% print(gcf,'-dpdf',[str,'.pdf'],'-r100');

