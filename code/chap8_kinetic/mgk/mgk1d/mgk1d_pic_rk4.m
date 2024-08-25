% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-08-31 16:25
% Delta-f PIC for ITG, electron adiabatic, with FLR, 1D, Dong 1992
% 17:19 can run, but not exact
% 16-09-03 LiYY modified version
% 16-09-04 09:10 seems can run, but not accurate, maybe normalization
% 16-09-20 09:42 recheck the normalization
% Rewrite 2017-02-06 23:06
% 17-02-10 17:02
% 17-02-14 08:53 try RK2->RK4, not yet
% 09:20 x=rand -> linspace, s
% 10:10 use Egt4(it)=mean(phig) instead of Egt1(it)=dphig(floor(ng/2)), seems ok
% 17-02-16 15:38 retry RK2->RK4

close all; clear; clc;
runtime=cputime;
irun=0;
nt=8000;
if(irun==0)
    dat=[0.45    -0.615095+0.263450i;
        0.50    -0.619550+0.245009i;
        0.55    -0.623357+0.222611i]; % From HanMK, Dong92
    id=1;
    kt=dat(id,1)/sqrt(2); % k_theta*rho_i
    wr=real(dat(id,2)); wi=imag(dat(id,2));
    s=1.0; q=1.0; tau=1.0; epsn=0.25; etai=2.5;
    tk=0.0; % theta_k

    wt0=(wr+1i*wi);
    wt=wt0*kt/epsn;

    np=40000*4;
    dt=0.0025;
    ng=128/1+0;

    xmax=4.0*pi; xmin=-xmax;
    L=xmax-xmin;
    dx=L/(ng-0);
    vt=1.0;
    vmax=6.0;

    zp=zeros(np,5); % x, vpar, vperp, p, w
    zp(:,1)=xmin+rand(np,1)*L;
%     deps=L/(np-1)/2;
%     zp(:,1)=linspace(xmin+deps,xmax-deps,np); % 17-02-14 08:55
%     zp(:,1)=xmin:L/(np-1):xmax;
    zp(:,2)=randn(np,1)*vt;
    zp(:,3)=sqrt(randn(np,1).^2+randn(np,1).^2)*vt;
    zp(:,4)=1.0;
    zp(:,5)=0.0001*exp(-zp(:,1).^2/1.5^2);
    czp=zp;

    pgmat=pginterp(zp,ng,dx,xmin,np);
    rhozero=sum(pgmat'*zp(:,4))/ng; % rhozero=np/ng;

    xg=xmin+(0:ng-1)'*dx;
    kyg2=kt^2.*(1+s^2*(xg-tk).^2);
    Gamma0=besseli(0,kyg2,1);
    wdi=-kt;
    wsi=-kt/epsn;
    qR_inv=1/q;
    hcoef=1./(1+1/tau-Gamma0);
    nta=1; ntb=nt;
else
    nta=ntb+1; ntb=nta+nt;
end

fd=2; % option of finite difference for dphi/dx
for it=nta:ntb
    
    zp(:,1)=mod(zp(:,1)-xmin+10*L,L)+xmin; % keep particles in [xmin,xmax]
    
    % The below 4 RK steps can be simplified, however, I just copy.
    % RK-4, 1st step
    ky=sqrt(kt^2.*(1+s^2*(zp(:,1)-tk).^2)); % kperp
    J0=besselj(0,ky.*zp(:,3));
    pgmat=pginterp(zp,ng,dx,xmin,np);
    rho=pgmat'*(zp(:,5).*J0)/rhozero;
    
    phig=rho.*hcoef; phig=smooth(phig);
    if(fd==1)
        phitmp=[phig(end);phig;phig(1)];
        dphig=(phitmp(3:end)-phitmp(1:ng))/dx/2;
    else
        phitmp=[phig(end-1);phig(end);phig;phig(1);phig(2)];
        dphig=(-phitmp(5:end)+phitmp(1:ng)+8*phitmp(4:end-1)-8*phitmp(2:ng+1))/dx/12;
    end
    dJ0=-besselj(1,ky.*zp(:,3)).*zp(:,3).*(zp(:,1)-tk)*s^2*kt^2./ky;
    phip=(pgmat*phig).*J0;
    dphip=(pgmat*dphig).*J0+(pgmat*phig).*dJ0;
    wD=wdi*(cos(zp(:,1))+s*(zp(:,1)-tk).*sin(zp(:,1))).*(zp(:,2).^2+zp(:,3).^2/2);
    wT=wsi*(1+(0.5*(zp(:,2).^2+zp(:,3).^2)-1.5)*etai);
    
    wdot1=-1i*wD.*zp(:,5)-qR_inv.*zp(:,2).*dphip+1i*(wT-wD).*phip;
    xdot1=qR_inv*zp(:,2);
    czp(:,1)=zp(:,1)+xdot1*dt*0.5;
    czp(:,5)=zp(:,5)+wdot1*dt*0.5;
    
    czp(:,1)=mod(czp(:,1)-xmin+10*L,L)+xmin;
    % RK-4, 2nd step
    ky=sqrt(kt^2.*(1+s^2*(czp(:,1)-tk).^2)); % kperp
    J0=besselj(0,ky.*czp(:,3));
    pgmat=pginterp(czp,ng,dx,xmin,np);
    rho=pgmat'*(czp(:,5).*J0)/rhozero;
    
    phig=rho.*hcoef; phig=smooth(phig);
    if(fd==1)
        phitmp=[phig(end);phig;phig(1)];
        dphig=(phitmp(3:end)-phitmp(1:ng))/dx/2;
    else
        phitmp=[phig(end-1);phig(end);phig;phig(1);phig(2)];
        dphig=(-phitmp(5:end)+phitmp(1:ng)+8*phitmp(4:end-1)-8*phitmp(2:ng+1))/dx/12;
    end
    dJ0=-besselj(1,ky.*czp(:,3)).*czp(:,3).*(czp(:,1)-tk)*s^2*kt^2./ky;
    phip=(pgmat*phig).*J0;
    dphip=(pgmat*dphig).*J0+(pgmat*phig).*dJ0;
    wD=wdi*(cos(czp(:,1))+s*(czp(:,1)-tk).*sin(czp(:,1))).*(czp(:,2).^2+czp(:,3).^2/2);
    wT=wsi*(1+(0.5*(czp(:,2).^2+czp(:,3).^2)-1.5)*etai);
    
    wdot2=-1i*wD.*czp(:,5)-qR_inv.*czp(:,2).*dphip+1i*(wT-wD).*phip;
    xdot2=qR_inv*czp(:,2);
    czp(:,1)=zp(:,1)+xdot2*dt*0.5;
    czp(:,5)=zp(:,5)+wdot2*dt*0.5;
    
    czp(:,1)=mod(czp(:,1)-xmin+10*L,L)+xmin;
    % RK-4, 3rd step
    ky=sqrt(kt^2.*(1+s^2*(czp(:,1)-tk).^2)); % kperp
    J0=besselj(0,ky.*czp(:,3));
    pgmat=pginterp(czp,ng,dx,xmin,np);
    rho=pgmat'*(czp(:,5).*J0)/rhozero;
    
    phig=rho.*hcoef; phig=smooth(phig);
    if(fd==1)
        phitmp=[phig(end);phig;phig(1)];
        dphig=(phitmp(3:end)-phitmp(1:ng))/dx/2;
    else
        phitmp=[phig(end-1);phig(end);phig;phig(1);phig(2)];
        dphig=(-phitmp(5:end)+phitmp(1:ng)+8*phitmp(4:end-1)-8*phitmp(2:ng+1))/dx/12;
    end
    dJ0=-besselj(1,ky.*czp(:,3)).*czp(:,3).*(czp(:,1)-tk)*s^2*kt^2./ky;
    phip=(pgmat*phig).*J0;
    dphip=(pgmat*dphig).*J0+(pgmat*phig).*dJ0;
    wD=wdi*(cos(czp(:,1))+s*(czp(:,1)-tk).*sin(czp(:,1))).*(czp(:,2).^2+czp(:,3).^2/2);
    wT=wsi*(1+(0.5*(czp(:,2).^2+czp(:,3).^2)-1.5)*etai);
    
    wdot3=-1i*wD.*czp(:,5)-qR_inv.*czp(:,2).*dphip+1i*(wT-wD).*phip;
    xdot3=qR_inv*czp(:,2);
    czp(:,1)=zp(:,1)+xdot3*dt;
    czp(:,5)=zp(:,5)+wdot3*dt;
    
    czp(:,1)=mod(czp(:,1)-xmin+10*L,L)+xmin;
    % RK-4, 4th step
    ky=sqrt(kt^2.*(1+s^2*(czp(:,1)-tk).^2)); % kperp
    J0=besselj(0,ky.*czp(:,3));
    pgmat=pginterp(czp,ng,dx,xmin,np);
    rho=pgmat'*(czp(:,5).*J0)/rhozero;
    
    phig=rho.*hcoef; phig=smooth(phig);
    if(fd==1)
        phitmp=[phig(end);phig;phig(1)];
        dphig=(phitmp(3:end)-phitmp(1:ng))/dx/2;
    else
        phitmp=[phig(end-1);phig(end);phig;phig(1);phig(2)];
        dphig=(-phitmp(5:end)+phitmp(1:ng)+8*phitmp(4:end-1)-8*phitmp(2:ng+1))/dx/12;
    end
    dJ0=-besselj(1,ky.*czp(:,3)).*czp(:,3).*(czp(:,1)-tk)*s^2*kt^2./ky;
    phip=(pgmat*phig).*J0;
    dphip=(pgmat*dphig).*J0+(pgmat*phig).*dJ0;
    wD=wdi*(cos(czp(:,1))+s*(czp(:,1)-tk).*sin(czp(:,1))).*(czp(:,2).^2+czp(:,3).^2/2);
    wT=wsi*(1+(0.5*(czp(:,2).^2+czp(:,3).^2)-1.5)*etai);
    
    wdot4=-1i*wD.*czp(:,5)-qR_inv.*czp(:,2).*dphip+1i*(wT-wD).*phip;
    xdot4=qR_inv*czp(:,2);
    zp(:,1)=zp(:,1)+(xdot1+2*xdot2+2*xdot3+xdot4)*dt/6.0;
    zp(:,5)=zp(:,5)+(wdot1+2*wdot2+2*wdot3+wdot4)*dt/6.0;
    
       
    EEk(it)=0.5*sum(zp(:,2).^2)/np*L; % kinetic energy
    EEf(it)=0.5*sum(dphig.^2)*dx; % potential energy
    Egt1(it)=dphig(floor(ng/2));
    Egt2(it)=mean(dphig);
    Egt3(it)=phig(floor(ng/2));
    Egt4(it)=mean(phig);
    
    %diagnosis
    if(it>=2 && mod(it,5)==0)        
        subplot(221);
        plot(xg,real(phig),xg,imag(phig),'LineWidth',2);
        title(['\rho, it=',num2str(it),'/',num2str(ntb)]);
        
        subplot(221); plot(2:it,real(Egt3(2:it)),2:it,...
            imag(Egt3(2:it)),'r--','LineWidth',2);
        title(['phi, it=',num2str(it),'/',num2str(ntb)]);
        subplot(222); plot(xg,real(phig),xg,imag(phig),'Linewidth',2);
        subplot(223); plot(2:it,real(Egt4(2:it)),2:it,...
            imag(Egt4(2:it)),'r--','LineWidth',2);
        subplot(224); plot(2:it,log(real(Egt4(2:it))),2:it,...
            log(imag(Egt4(2:it))),'r--','LineWidth',2);
        drawnow;
    end
end
runtime=cputime-runtime;

%%
close all;
h = figure('Unit','Normalized','position',...
    [0.02 0.1 0.6 0.7],'DefaultAxesFontSize',15);
Egt=Egt4;
t=(1:ntb)*dt;
subplot(221);
plot(t,real(Egt),t,imag(Egt),'--','LineWidth',2);
title(['k_\theta=',num2str(kt),', s=',num2str(s),', q=',num2str(q),...
    ', \epsilon_n=',num2str(epsn),', \eta_i=',num2str(etai),...
    ', \tau=',num2str(tau)]);
xlabel(['t, runtime=',num2str(runtime),'s']); 
subplot(222);
% Find the corresponding indexes of the extreme max values
% lndE=log(sqrt(real((EEf(1:nt))))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
lndE=log(real((Egt(1:ntb))));
lndEi=log(imag((Egt(1:ntb))));
lndEa=log(abs((Egt(1:ntb))));
it0=floor(ntb*3/20); it1=floor(ntb*19.8/20);
yy=lndE(it0:it1);
plot(t,lndE,'LineWidth',2); hold on;
pltwr=2; %
yy1=smooth(real(yy),15); % 16-12-12 10:26
if(pltwr==1)
    extrMaxIndex = find(diff(sign(diff(yy)))==-2)+1;
    t1=t(it0+extrMaxIndex(1));t2=t(it0+extrMaxIndex(end));
    y1=yy(extrMaxIndex(1));y2=yy(extrMaxIndex(end));
    plot([t1,t2],[y1,y2],'r*--','LineWidth',2);
    nw=length(extrMaxIndex)-1;
    omega=pi/((t2-t1)/nw);
elseif(pltwr==2)
    [pks,locs]=findpeaks(yy1,'minpeakdistance',200);
    t1=t(it0+locs(1));t2=t(it0+locs(end));
    y1=yy(locs(1));y2=yy(locs(end));
    nw=length(locs)-1;
    omega=pi/((t2-t1)/nw);
else
    t1=t(it0);t2=t(it1);y1=yy(1);y2=yy(end);
    nw=0; omega=0;
end
plot([t1,t2],[y1,y2],'r*--','LineWidth',2);
for jp=1:nw+1
    hold on; plot(t(it0+locs(jp)),yy(locs(jp)),'g+');
end
gammas=(real(y2)-real(y1))/(t2-t1);

title(['(b) \omega^S=',num2str(omega+1i*gammas),...
    ', nw=',num2str(nw)]);
axis tight; xlabel(['\omega^T=',num2str(wt,3)]);

subplot(223);
% plot(xg,abs(phig),'LineWidth',2);
phig=phig/phig(floor(ng/2));
plot(xg,real(phig),xg,imag(phig),'--','LineWidth',2);

th_pic=xg;
phi_pic=phig;

if(id==1)
    run phi_hd7_dat; hold on;
    plot(th_hd7,real(phi_hd7),'k:',...
        th_hd7,imag(phi_hd7),'m:','Linewidth',2);
end
subplot(224);
scatter(zp(:,2),zp(:,3),abs(zp(:,5)),abs(zp(:,5)));

print(gcf,'-dpng',['mgk1d_pic_rk4_s=',num2str(s),',q=',num2str(q),...
    ',tau=',num2str(tau),',eta=',num2str(etai),',epsn=',num2str(epsn),...
    ',kt=',num2str(kt),',ng=',num2str(ng),...
    ',np=',num2str(np),',dt=',num2str(dt),...
    ',xmax=',num2str(xmax),',fd=',num2str(fd),'.png']);
