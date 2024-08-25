% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-04-20 13:34
% Delta-f PIC for IAW, electron adiabatic
% Test OK, 2016-04-20 20:43
close all; 
clear; clc;
dat2 = [1,2.04590486664825,-0.851330458171737; % sqrt(2)
    2,2.37997049716181,-0.568627785772305;
    3,2.62647608379667,-0.411054266743725;
    4,2.83132357729409,-0.306718776560365;
    5,3.01110302739500,-0.232235873041931;
    6,3.17397204092200,-0.176905646874478;
    7,3.32463253189007,-0.134880072940063;
    8,3.46607414950534,-0.102578518320425;
    9,3.60032687530566,-0.0776243048338077;
    10,3.72883249481195,-0.0583410544450318];
runtime=cputime;
k=1.0; nfilter=1; deltaf=1; linear=0;
id=2;
tau=dat2(id,1); wr=dat2(id,2); wi=dat2(id,3);

np = 160000/4;
dt=0.01/2;
nt=1000;
ng=16;
L=2*pi/k;
dx=L/ng;
% vt=1.0; % vt=sqrt(T/m)
vmax=6.0;

zp=zeros(np,4); % x, v, p, w
zp(:,1)=rand(np,1)*L;
zp(:,2)=(rand(np,1)*2-1)*vmax;
zp(:,3)=exp(-zp(:,2).^2/2)/sqrt(2*pi);
zp(:,4)=0.00001*zp(:,3).*cos(k*zp(:,1));
pgmat=pginterp(zp,ng,dx,np);
rhozero=sum(pgmat'*zp(:,3))/ng; % rhozero=np/ng; 

xg=(0:ng-1)'*dx;

%%
figure;
for it=1:nt
    
    zp(:,1)=mod(zp(:,1)+L,L); % keep particles in [0,L]
    
    % RK-2, 1st step
    pgmat=pginterp(zp,ng,dx,np);
    rho=pgmat'*zp(:,4)/rhozero;
    
    if(nfilter==1)
        rhof=fft(rho)/ng;
        af1=real(rhof(2)+rhof(ng)); bf1=imag(-rhof(2)+rhof(ng));
        Eg=-tau*(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg));
    else
        phi=tau*[rho(end);rho;rho(1)];
        Eg=-(phi(3:end)-phi(1:ng))/dx/2;
    end
    
    vdot=pgmat*Eg;
    wdot=(zp(:,3)-(1-linear)*zp(:,4)).*vdot.*zp(:,2);
    czp(:,1)=zp(:,1)+zp(:,2)*dt/2;
    czp(:,2)=zp(:,2)+vdot*dt/2;
    czp(:,3)=zp(:,3);
    czp(:,4)=zp(:,4)+wdot*dt/2;
    
    czp(:,1)=mod(czp(:,1)+L,L);
    % RK-2, 2nd step
    pgmat=pginterp(czp,ng,dx,np);
    rho=pgmat'*czp(:,4)/rhozero;
    if(nfilter==1)
        rhof=fft(rho)/ng;
        af1=real(rhof(2)+rhof(ng)); bf1=imag(-rhof(2)+rhof(ng));
        Eg=-tau*(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg));
    else
        phi=tau*[rho(end);rho;rho(1)];
        Eg=-(phi(3:end)-phi(1:ng))/dx/2;
    end
    
    vdot=pgmat*Eg;
    wdot=(czp(:,3)-(1-linear)*czp(:,4)).*vdot.*czp(:,2);
    zp(:,1)=zp(:,1)+czp(:,2)*dt;
    zp(:,2)=zp(:,2)+vdot*dt;
    zp(:,3)=zp(:,3);
    zp(:,4)=zp(:,4)+wdot*dt;
       
    %diagnosis
    subplot(121);
    plot(xg,rho,'LineWidth',2);
    title(['xhs, \rho, it=',num2str(it),'/',num2str(nt)]);
    subplot(122);
    plot(xg,Eg,'LineWidth',2);
    title('Eg');
    
    EEk(it)=0.5*sum(zp(:,2).^2)/np*L; % kinetic energy
    EEf(it)=0.5*sum(Eg.^2)*dx; % potential energy
    drawnow;
end

runtime=cputime-runtime;

%%
close all;
h = figure('Unit','Normalized','position',...
    [0.02 0.2 0.6 0.4],'DefaultAxesFontSize',15);
t=(1:nt)*dt;
EEf1=EEf/tau;
subplot(121);
plot(t,EEf1,t,EEk,t,EEk+EEf1,'r--','LineWidth',2);
title(['(a) \tau=',num2str(tau),', runtime=',num2str(runtime)]);
xlabel('t');
subplot(122);
% Find the corresponding indexes of the extreme max values
lndE=log(sqrt(real((EEf(1:nt))))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
% lndE=log(real((Egt(1:nt))));
it0=floor(nt*1/20); it1=floor(nt*18/20);
yy=lndE(it0:it1);
plot(t,lndE,'LineWidth',2);
hold on;
extrMaxIndex = find(diff(sign(diff(yy)))==-2)+1;
t1=t(it0+extrMaxIndex(1));t2=t(it0+extrMaxIndex(end));
y1=yy(extrMaxIndex(1));y2=yy(extrMaxIndex(end));
plot([t1,t2],[y1,y2],'r*--','LineWidth',2);
nw=length(extrMaxIndex)-1;
omega=pi/((t2-t1)/nw);
gammas=(real(y2)-real(y1))/(t2-t1);
title(['(b) \omega^S=',num2str(omega),', \gamma^S=',num2str(gammas),...
    ', nw=',num2str(nw)]);
axis tight;
xlabel(['\omega^T=',num2str(wr+1i*wi,3)]);

str=['pic1d_df_iaw_k=',num2str(k),',ng=',num2str(ng),'np=',num2str(np),',nt=',num2str(nt),...
    ',dt=',num2str(dt),',tau=',num2str(tau),',nfilter=',num2str(nfilter)];
print(gcf,'-dpng',[str,'_history.png']);
