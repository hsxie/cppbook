% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-04-21 10:08
% Delta-f PIC for bump-on-tail
% ok for Landau damping, but not ok for beam, why? 16-04-21 16:34
% beam ok, because nt*dt and vd too small at previous test 16-04-22 15:52
%
irun=0;

nt=1000;
runtimea=cputime;
close all; clc; 
if(irun==0)
    clear EEk  EEf Egt t;

    data=[0.0500000000000000,0.219444454284370,0.0322421790271008;
    0.100000000000000,0.425776264847492,0.0708336313405495;
    0.150000000000000,0.603780584878508,0.111524902450115;
    0.200000000000000,0.750156727061233,0.138177445176053;
    0.250000000000000,0.875204244169648,0.142060277851812;
    0.300000000000000,0.987304619094534,0.124378064165398;
    0.350000000000000,1.08989854927862,0.0885726360304440;
    0.400000000000000,1.18369937208857,0.0379333673216506]; % nb=0.1, vd=5.0
    id=6;
    k=data(id,1);
    wr=data(id,2); wi=data(id,3);

    runtime=0;
    nfilter=0; deltaf=1; nonlinear=1;

    np = 5000;
%     nt=2000;
    dt=0.05;
    ng=32*1;
    L=2*pi/k;
    dx=L/ng;
    % vt=1.0; % vt=sqrt(T/m)
    vmax=8.0;

    nb=0.1; vdb=5.0;
    zp=zeros(np,4); % x, v, p, w
    zp(:,1)=rand(np,1)*L;
    zp(:,2)=(rand(np,1)*2-1)*vmax+0*vdb/2;
    elec=2;
    if(elec==1)
        zp(:,3)=exp(-zp(:,2).^2/2)/sqrt(2*pi);
    else
        zp(:,3)=(1-nb)*exp(-zp(:,2).^2/2)/sqrt(2*pi)+...
            nb*exp(-(zp(:,2)-vdb).^2/2)/sqrt(2*pi);
    end
    zp(:,4)=0.00001*zp(:,3).*sin(k*zp(:,1));
%     zp(:,4)=0.0001*(rand(np,1)-0.5);
    rhozero=sum(zp(:,3))/ng; % rhozero=np/ng;

    xg=(0:ng-1)'*dx;
    nta=1; ntb=nt;
else
    nta=it+1; ntb=nta+nt-1;
end

%%
figure('Unit','Normalized','position',...
    [0.02 0.1 0.65 0.65],'DefaultAxesFontSize',15);
for it=nta:ntb
    
    zp(:,1)=mod(zp(:,1)+10*L,L); % keep particles in [0,L]
    
    % RK-2, 1st step
    pgmat=pginterp(zp,ng,dx,np);
    rho=pgmat'*zp(:,4)/rhozero;
    
    if(nfilter==1)
        rhof=fft(rho)/ng;
        af1=real(rhof(2)+rhof(ng)); bf1=imag(-rhof(2)+rhof(ng));
        Eg=-(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg))/k;
    else
        Eg=cumsum(rho)*dx;
        Eg=Eg-mean(Eg);
    end
    
    vdot=pgmat*Eg;
    if(elec==1)
        wdot=(zp(:,3)-nonlinear*zp(:,4)).*vdot.*zp(:,2);
    else
        dlnf0dv=-(zp(:,2).*(1-nb).*exp(-zp(:,2).^2/2)+(zp(:,2)-...
            vdb).*nb.*exp(-(zp(:,2)-vdb).^2/2))./((1-nb).*exp(-zp(:,2).^2/2)+...
            nb.*exp(-(zp(:,2)-vdb).^2/2));
        wdot=-(zp(:,3)-nonlinear*zp(:,4)).*vdot.*dlnf0dv;
    end
    czp(:,1)=zp(:,1)+zp(:,2)*dt/2;
    czp(:,2)=zp(:,2)+nonlinear*vdot*dt/2;
    czp(:,3)=zp(:,3);
    czp(:,4)=zp(:,4)+wdot*dt/2;
    
    czp(:,1)=mod(czp(:,1)+10*L,L);
    % RK-2, 2nd step
    pgmat=pginterp(czp,ng,dx,np);
    rho=pgmat'*czp(:,4)/rhozero;
    if(nfilter==1)
        rhof=fft(rho)/ng;
        af1=real(rhof(2)+rhof(ng)); bf1=imag(-rhof(2)+rhof(ng));
        Eg=-(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg))/k;
    else
        Eg=cumsum(rho)*dx;
        Eg=Eg-mean(Eg);
    end
    
    vdot=pgmat*Eg;
    if(elec==1)
        wdot=(czp(:,3)-nonlinear*czp(:,4)).*vdot.*czp(:,2);
    else
        dlnf0dv=-(czp(:,2).*(1-nb).*exp(-czp(:,2).^2/2)+(czp(:,2)-...
            vdb).*nb.*exp(-(czp(:,2)-vdb).^2/2))./((1-nb).*exp(-czp(:,2).^2/2)+...
            nb.*exp(-(czp(:,2)-vdb).^2/2));
        wdot=-(czp(:,3)-nonlinear*czp(:,4)).*vdot.*dlnf0dv;
    end
    zp(:,1)=zp(:,1)+czp(:,2)*dt;
    zp(:,2)=zp(:,2)+nonlinear*vdot*dt;
    zp(:,3)=zp(:,3);
    zp(:,4)=zp(:,4)+wdot*dt;

    EEk(it)=0.5*sum(zp(:,2).^2)/np*L; % kinetic energy
    EEf(it)=0.5*sum(Eg.^2)*dx; % potential energy
    Egt(it)=Eg(floor(ng/2));
    
    % diagnosis
    subplot(231);
    plot(xg,rho,'LineWidth',2);
    title(['\rho, it=',num2str(it),'/',num2str(ntb)]);
    xlim([min(xg),max(xg)]);
    subplot(232);
    plot(xg,Eg,'LineWidth',2);
    title('Eg');
    xlim([min(xg),max(xg)]);
    subplot(233);
    plot(zp(:,2),zp(:,3),'.'); xlabel('v'); ylabel('p');
    xlim([-1.1*vmax,1.1*vmax]);
    subplot(234);
    scatter(zp(:,1),zp(:,2),zp(:,3),zp(:,3)+zp(:,4));
    xlabel('x'); ylabel('v');
    ylim([-1.1*vmax,1.1*vmax]); xlim([min(xg),max(xg)]);
    subplot(235);
    plot((1:it)*dt,real(log(EEf)),'Linewidth',2); 
    xlabel('t'); ylabel('log(Ef)');
    xlim([0,it*dt]);
    subplot(236);
    scatter(zp(:,1),zp(:,2),zp(:,3),zp(:,4));
    ylim([-1.1*vmax,1.1*vmax]); xlim([min(xg),max(xg)]);
    xlabel('x'); ylabel('v');
    
    drawnow;
end
runtime=runtime+cputime-runtimea;

str=['pic1d_df_beam_k=',num2str(k),',ng=',num2str(ng),'np=',...
    num2str(np),',nt=',num2str(ntb),',dt=',num2str(dt),...
    ',nl=',num2str(nonlinear),',nfilter=',num2str(nfilter)];

print(gcf,'-dpng',[str,'_snap.png']);
%%
% close all;
h = figure('Unit','Normalized','position',...
    [0.02 0.2 0.6 0.4],'DefaultAxesFontSize',15);
t=(1:ntb)*dt;
EEf1=EEf;
subplot(121);
plot(t,EEf1,t,EEk,t,EEk+EEf1,'r--','LineWidth',2);
title(['(a) k=',num2str(k),', runtime=',num2str(runtime)]);
xlabel('t');

subplot(122);
% Find the corresponding indexes of the extreme max values
% lndE=log(sqrt(real((EEf(1:nt))))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
lndE=log(real((Egt(1:ntb))));
it0=floor(ntb*7/20); it1=floor(ntb*17/20);
yy=real(lndE(it0:it1));
plot(t,lndE,'LineWidth',2); hold on;
[pks,locs] = findpeaks(yy);
t1=t(it0+locs(1));t2=t(it0+locs(end));
y1=yy(locs(1));y2=yy(locs(end));
plot([t1,t2],[y1,y2],'r*--','LineWidth',2);
nw=length(locs)-1;
omega=pi/((t2-t1)/nw);
gammas=(real(y2)-real(y1))/(t2-t1);
title(['(b) \omega^S=',num2str(omega),', \gamma^S=',num2str(gammas),...
    ', nw=',num2str(nw)]);
axis tight;
xlabel(['\omega^T=',num2str(wr+1i*wi,3)]);

str=['pic1d_df_beam_k=',num2str(k),',ng=',num2str(ng),'np=',...
    num2str(np),',nt=',num2str(ntb),',dt=',num2str(dt),...
    ',nl=',num2str(nonlinear),',nfilter=',num2str(nfilter)];
print(gcf,'-dpng',[str,'_history.png']);
