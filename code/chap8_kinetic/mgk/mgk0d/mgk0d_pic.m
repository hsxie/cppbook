% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-09-01 00:02
% Delta-f PIC for entropy, ion & electron both kinetic, with FLR, 0D
% 16-10-10 20:24-21:55 rewrite, can run and agree with D.R.

% to add kpara

close all; clear; clc;

runtimea=cputime;
irun=0;
nt=500;
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
    Gamma0i=besseli(0,ki^2,1);
    Gamma0e=besseli(0,ke^2,1);
    G0coef=1.0/((1.0-Gamma0i)+(1.0-Gamma0e)/tau);

    ms=[1,1/mi];
    qs=[1,-1];
    ks=[ki,ke]*sqrt(1);
    wds=[wdi,wde];
    Ti=1.0; Te=Ti*tau;
    taus=[1,tau];
    Ts=[Ti,Te];
    vts=sqrt(Ts./ms);

    np = 40000*1/1; dt=0.02;
    vmax=6.0;

    ns=2;
    zp=zeros(ns,np,3); % vpar, vperp, w
    for is=1:ns
        zp(is,:,1)=randn(np,1);
        zp(is,:,2)=sqrt(randn(np,1).^2+randn(np,1).^2);
        zp(is,:,3)=qs(is)*0.00001;
%         zp(is,:,3)=0.00001;
        rhozero(is)=np; % rhozero=np/ng; 
    end
    czp=zp;

    tt=[];
    phit=[];
    nta=1; ntb=nt;
else
    nta=ntb+1; ntb=nta+nt;
end

%%
figure;
for it=nta:ntb
    
    % RK-2, 1st step
    rho=0;
    for is=1:ns
        rho=rho+sum(zp(is,:,3))/rhozero(is)/taus(is);
    end
    phi=G0coef*rho;
    
    for is=1:ns
        wdot=-1i*(wds(is)*(zp(is,:,1).^2+0.5*zp(is,:,2).^2)).*zp(is,:,3)...
            -1i*phi*(wds(is)*(zp(is,:,1).^2+0.5*zp(is,:,2).^2)...
            -wds(is)*(kapn+(0.5*(zp(is,:,1).^2+zp(is,:,2).^2)-1.5)*kapt)...
            ).*besselj(0,ks(is)*zp(is,:,2)).^2;
        czp(is,:,3)=zp(is,:,3)+wdot*dt/2;
    end
    
    % RK-2, 2nd step
    rho=0;
    for is=1:ns
        rho=rho+sum(czp(is,:,3))/rhozero(is)/taus(is);
    end
    phi=G0coef*rho;
    
    for is=1:ns
        wdot=-1i*(wds(is)*(czp(is,:,1).^2+0.5*czp(is,:,2).^2)).*czp(is,:,3)...
            -1i*phi*(wds(is)*(czp(is,:,1).^2+0.5*czp(is,:,2).^2)...
            -wds(is)*(kapn+(0.5*(czp(is,:,1).^2+czp(is,:,2).^2)-1.5)*kapt)...
            ).*besselj(0,ks(is)*czp(is,:,2)).^2;
        zp(is,:,3)=zp(is,:,3)+wdot*dt;
    end
    
    tt=[tt,it*dt];
    phit=[phit,phi];
    %diagnosis
    if(mod(it,100)==1)
        subplot(121);
        plot(tt,real(phit),tt,imag(phit),'r--','LineWidth',2);
        title(['\phi, it=',num2str(it),'/',num2str(ntb)]);

        subplot(122);
        plot(tt,log(real((phit))),'LineWidth',2); box on;
        xlabel('t'); ylabel('log(phit)');
    
        drawnow;
    end
end

runtime=runtime+cputime-runtimea;

%%
close all;
h = figure('Unit','Normalized','position',...
    [0.02 0.2 0.6 0.4],'DefaultAxesFontSize',15);
ntc=length(phit);
t=(1:ntc)*dt;
subplot(121);
plot(tt,real(phit),tt,imag(phit),'r--','LineWidth',2);
title(['(a) k_\perp=',num2str(ky),...
    ', \epsilon_n=',num2str(epsn),', \eta=',num2str(eta),...
    ', \tau=',num2str(tau)]);
xlabel(['t',', runtime=',num2str(runtime,4),'s']); 
subplot(122);
% Find the corresponding indexes of the extreme max values
% lndE=log(sqrt(real((EEf(1:nt))))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
lndE=log(real((phit(1:ntc))));
lndEi=log(imag((phit(1:ntc))));
lndEa=log(abs((phit(1:ntc))));
it0=floor(ntc*5.5/20); it1=floor(ntc*19.9/20);
yy=lndE(it0:it1);
plot(t,lndE,t,lndEi,'g--',t,lndEa,'k:','LineWidth',2);
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
axis tight;
xlabel(['\omega^T=',num2str(wr+1i*wi,4)]);

str=['gkpic0d_entropy_ky=',num2str(ky), ',etai=',num2str(eta),...
     ',epsn=',num2str(epsn),',tau=',num2str(tau),...
     ',np=',num2str(np),',nt=',num2str(nt),',dt=',num2str(dt)];
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 4.0]);
print(gcf,'-dpng',[str,'_history_2.png']);
