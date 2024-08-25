% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-08-30 10:24
% Delta-f PIC for ITG, electron adiabatic, with FLR, 1D
% 20:06 can run, but nring not convergent
% 16-10-01 00:02 modify to kapn, kapt, seems can run but not exact to DR
% 16-10-02 17:55 accurate for DR now, Cummings' nring=4 accurate for kperp<1.0 

close all; clear; clc;

irun=0;

nt=200*1;
runtimea=cputime;

if(irun==0)
    
%     dat=[0.1000    0.0822    0.0716;
%     0.2000    0.1676    0.1342;
%     0.3000    0.2575    0.1807;
%     0.4000    0.3503    0.2072;
%     0.5000    0.4428    0.2135;
%     0.6000    0.5312    0.2026;
%     0.7000    0.6125    0.1803;
%     0.8000    0.6847    0.1537;
%     0.9000    0.7463    0.1294;
%     1.0000    0.7954    0.1115;
%     1.1000    0.8312    0.0999;
%     1.2000    0.8547    0.0924]; % kapn=0.0, kapt=2.5
    dat=[0.1000    0.0219    0.0628;
    0.2000    0.0522    0.1208;
    0.3000    0.0959    0.1687;
    0.4000    0.1525    0.2019;
    0.5000    0.2179    0.2175;
    0.6000    0.2863    0.2159;
    0.7000    0.3524    0.2000;
    0.8000    0.4131    0.1745;
    0.9000    0.4673    0.1445;
    1.0000    0.5157    0.1148;
    1.1000    0.5587    0.0899;
    1.2000    0.5957    0.0723;
    2.0 0 0]; % kapn=1.0, kapt=2.5
    id=10;
    k=dat(id,1); wr=dat(id,2); wi=dat(id,3); 
    
    nfilter=1; deltaf=1; nonlinear=0;
    
    runtime=cputime;
    
%     kapn=0.0; kapt=2.5;
    kapn=1.0; kapt=2.5;
    tau=1.0;
    epsn=0.2; 
    wd=2*epsn;
%     wd=0;
    
    Gamma0=besseli(0,k^2,1); % Gamma0=1.0; %

    np = 4000*1/1;
    nring=8*1;
    np2=np*nring;
    dt=0.2/k;
    ng=16*4;
    L=2*pi/k;
    dx=L/ng;
    T=0.5; % 
    m=1;
    vt=sqrt(T/m);
    vmax=6.0;
    

    zp=zeros(np,5); % x, vpar, vperp, p, w
    zpring=zeros(np2,1);
    zp(:,1)=rand(np,1)*L;
    zp(:,2)=randn(np,1)*vt;
    zp(:,3)=sqrt(randn(np,1).^2+randn(np,1).^2)*vt;
    zp(:,4)=1.0;
    % zp(:,5)=0.00001*zp(:,3).*cos(k*zp(:,1));
    zp(:,5)=0.00001*zp(:,4).*cos(k*zp(:,1));
    % pgmat=pginterp(zp,ng,dx,np);
    % rhozero=sum(pgmat'*zp(:,4))/ng; % rhozero=np/ng; 
    rhozero=np/ng; 

    xg=(0:ng-1)'*dx;
    tt=[];phit=[];phi2t=[];
    nta=1; ntb=nt;
else
    nta=it+1; ntb=nta+nt-1;
end

%%
figure('Unit','Normalized','position',...
    [0.02 0.1 0.65 0.65],'DefaultAxesFontSize',15);
for it=nta:ntb
    
    zp(:,1)=mod(zp(:,1)+10*L,L); % keep particles in [0,L]
    zpring = zctop(zp,np,nring,L);
    
    % RK-2, 1st step
    zcpmat = zcpinterp(zpring,ng,dx,np,nring,L);
    rho=zcpmat'*zpring(:,2)/rhozero;
    
    if(nfilter==1)
        rhof=fft(rho)/ng;
        af1=real(rhof(2)+rhof(ng)); bf1=imag(-rhof(2)+rhof(ng));
%         Eg=-tau*(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg));
%         Eg=-k*(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg))/(1/tau+k^2);
        Eg=-k*(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg))/(1+tau-Gamma0);
%         Eg=k*(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg))/(1+tau-Gamma0);
    else
        phi=tau*[rho(end);rho;rho(1)];
        Eg=-(phi(3:end)-phi(1:ng))/dx/2;
    end
    
    Epring=zcpmat*Eg;
    Ep = Eptoc(Epring,np,nring);
    
    vdot=Ep;
    wdot=-vdot.*((kapn+(zp(:,2).^2+zp(:,3).^2-1.5)*kapt)- ...
        wd*(zp(:,2).^2+0.5*zp(:,3).^2));
    czp(:,1)=zp(:,1)+wd*(zp(:,2).^2+0.5*zp(:,3).^2)*dt/2;
    czp(:,5)=zp(:,5)+wdot*dt/2;
    
    czp(:,1)=mod(czp(:,1)+10*L,L);
    czpring = zctop(czp,np,nring,L);
    % RK-2, 2nd step    
    czcpmat = zcpinterp(czpring,ng,dx,np,nring,L);
%     rho=zcpmat'*czpring(:,2)/rhozero/nring;
    rho=zcpmat'*czpring(:,2)/rhozero;
    if(nfilter==1)
        rhof=fft(rho)/ng;
        af1=real(rhof(2)+rhof(ng)); bf1=imag(-rhof(2)+rhof(ng));
        Eg=-k*(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg))/(1+tau-Gamma0);
    else
        phi=tau*[rho(end);rho;rho(1)];
        Eg=-(phi(3:end)-phi(1:ng))/dx/2;
    end
    
    Epring=zcpmat*Eg;
    Ep = Eptoc(Epring,np,nring);
    
    vdot=Ep;
    wdot=-vdot.*((kapn+(zp(:,2).^2+zp(:,3).^2-1.5)*kapt)- ...
        wd*(zp(:,2).^2+0.5*zp(:,3).^2));
    zp(:,1)=zp(:,1)+wd*(zp(:,2).^2+0.5*zp(:,3).^2)*dt;
    zp(:,5)=zp(:,5)+wdot*dt;
       
    %diagnosis
    subplot(121);
    plot(xg,rho,'LineWidth',2);
    title(['xhs, \rho, it=',num2str(it),'/',num2str(nt)]);
    subplot(122);
    plot(xg,Eg,'LineWidth',2);
    title('Eg');
    
    EEk(it)=0.5*sum(zp(:,2).^2)/np*L; % kinetic energy
    EEf(it)=0.5*sum(Eg.^2)*dx; % potential energy
    Egt(it)=Eg(floor(ng/2));
    drawnow;
end

runtime=cputime-runtime;

%%
close all;
h = figure('Unit','Normalized','position',...
    [0.02 0.2 0.6 0.4],'DefaultAxesFontSize',15);
t=(1:ntb)*dt;
EEf1=EEf/tau;
subplot(121);
% plot(t,EEf1,t,EEk,t,EEk+EEf1,'r--','LineWidth',2);
plot(xg,Eg,'LineWidth',2);
title(['(a) \tau=',num2str(tau),', k=',num2str(k),...
    ', runtime=',num2str(runtime,3),'s']);
xlabel('t'); axis tight;
subplot(122);
% Find the corresponding indexes of the extreme max values
% lndE=log(sqrt(real((EEf(1:nt))))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
lndE=log(real((Egt(1:ntb))));
it0=floor(ntb*4/20); it1=floor(ntb*19.8/20);
yy=lndE(it0:it1);
plot(t,lndE,'LineWidth',2);
hold on;
extrMaxIndex = find(diff(sign(diff(yy)))==-2)+1;
t1=t(it0+extrMaxIndex(1));t2=t(it0+extrMaxIndex(end));
y1=yy(extrMaxIndex(1));y2=yy(extrMaxIndex(end));
% t1=t(it0);t2=t(it1);
% y1=yy(1);y2=yy(end);
plot([t1,t2],[y1,y2],'r*--','LineWidth',2);
nw=length(extrMaxIndex)-1;
omega=pi/((t2-t1)/nw);
% omega=0; nw=1;
gammas=(real(y2)-real(y1))/(t2-t1);
title(['(b) \omega^S=',num2str(omega),', \gamma^S=',num2str(gammas),...
    ', nw=',num2str(nw)]);
axis tight;
xlabel(['\omega^T=',num2str(wr+1i*wi,3)]);

str=['gkpic1d_itg_flr_k=',num2str(k),',ng=',num2str(ng),...
    ',kapn=',num2str(kapn),',kapt=',num2str(kapt),',wd=',num2str(wd),...
    ',np=',num2str(np),',nring=',num2str(nring),',ntb=',num2str(ntb),...
    ',dt=',num2str(dt),',tau=',num2str(tau),',nfilter=',num2str(nfilter)];
print(gcf,'-dpng',[str,'_history.png']);
