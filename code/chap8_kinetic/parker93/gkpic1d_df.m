% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-07-02 15:28
% Gyro-kinetic Delta-f PIC 1D drift wave
% Ref: a. Parker&Lee1993; b. Lee lecture notes
% 1. Landau damping, only electron, ok, 16:00
% 2. IAW, both electron and ion, 16-07-08 09:07, can run, but not exact yet
% (16-08-17 add: It's exact, but for electron branch due to high frequency)
% 3. Parker & Lee, 1993, PFB, 16-07-15 12:40
% 16-08-16 21:15, ok for kappan=0, agree with gk1d_dr.m, but not exact for
% kappan \neq 0
% 16-08-18 14:07 kappat \neq 0, kappen=0, agree with gk1d_dr_with_kt.m
% 16-09-18 00:01 Zhang Yi found the bug: RK-2 2nd step fft(rhos) should be
% fft(rho)!! Both linear and nonlinear similar to Parker93 now.
% 16-10-02 23:43 with also only ion, test nonlinear, not ok
% 16-10-03 12:26, cal transport coeffients, Ref Lee&Tang88, Eqs.22-25
% 16-10-24 10:56 add fobanacci loading

close all; clear; clc;
irun=0;

nt=5000;
runtimea=cputime;
close all; clc;
if(irun==0)
    
    clear EEk  EEf Egt t;
%     data = [0.5   0.4915 0.498 ;
%         1.0 0.675 0.482]; %  kappat=5.0;
%     data = [0.38   0.391 0.131 ;
%         1.0 0.611 -0.056]; %  kappat=1.0;
%     data = [0.5   0.291 0.176 ;
%         1.0 0.326 0.191]; %  kappan=1.0;
    data = [1.2   0.0 0.0;
        0.5   0.0804 0.0074;
        0.8   0.0966 0.0131;
        0.98   0.0988 0.0138;
        1.0 0.0989 0.0138]; %  kappan=0.2;

%     data = [0.5   -0.0146 0.0217 ;
%         1.0 0.0253 0.0362];
    
    id=3;
    k=data(id,1);
    wr=data(id,2); wi=data(id,3); Gamma0=besseli(0,k^2,1);
%     theta=0.01; kappan=0.0; kappat=0.2;
%     theta=0.01; kappan=0.0; kappat=1.0;
    theta=0.01; kappan=0.2; kappat=0.0;
%     theta=0.01; kappan=0.4; kappat=0.4;
    
    runtime=0;
    nfilter=1; deltaf=1; %nonlinear=1;
    nload=2;

    
    only_ion=2;
    if(only_ion==1)
        ns=1; % adiabatic electron 
        dt=1.0/1/4;
    else
        ns=2; % two species
        dt=0.2;
    end
%     ms=[1,1/25];
    ms=[1,1/1837];
%     ms=[1837,1];
    qs=[1,-1];
%     nonlinears=[0,0];
%     nonlinears=[1,1];
    nonlinears=[0,1]; % linear ion in Parker93
    np = 10000*4*1;
    ng=32*2/1;
    L=2*pi/k;
    dx=L/ng;
    tau=1.0;
    % vt=1.0; % vt=sqrt(T/m)
    Ts=[1.0,1.0];
%     Ts=[2,2];
    vts=sqrt(Ts./ms);
%     vts=sqrt(2*Ts./ms);
    vmax=5.0;

    if(nload==1)
        zp=zeros(ns,np,4); % x, v, p, w
        czp=zp;
        for is=1:ns
            zp(is,:,1)=rand(np,1)*L;
            %         zp(is,:,2)=randn(np,1);
            %         zp(is,:,2)=randn(np,1)*vts(is)/sqrt(2); % randn is exp(-v^2/2)
            zp(is,:,2)=randn(np,1)*vts(is);
            zp(is,:,3)=1.0;
            zp(is,:,4)=0.0001*zp(is,:,3).*sin(k*zp(is,:,1));
            rhozero(is)=sum(zp(is,:,3))/ng; % rhozero=np/ng;
        end
    else % fobanacci loading        
        M=22; % M=23;
        [xp,vp,np]=fun_fobanacci(M);
        zp=zeros(ns,np,4); % x, v, p, w
        czp=zp;
        for is=1:ns
            zp(is,:,1)=xp*L;
            zp(is,:,2)=vp*vts(is);
            zp(is,:,3)=1.0;
            zp(is,:,4)=0.0001*zp(is,:,3).*sin(k*zp(is,:,1));
            rhozero(is)=sum(zp(is,:,3))/ng; % rhozero=np/ng;
        end
    end

    xg=(0:ng-1)'*dx;
    nta=1; ntb=nt;
%     EEk=zeros(ntb,1); EEf=zeros(ntb,1); Egt=zeros(ntb,1);
%     denflux=zeros(ns,ntb); engflux=zeros(ns,ntb);
else
    nta=it+1; ntb=nta+nt-1;
end

%
figure('Unit','Normalized','position',...
    [0.02 0.06 0.65 0.8],'DefaultAxesFontSize',15);
for it=nta:ntb
    
    zp(:,:,1)=mod(zp(:,:,1)+10*L,L); % keep particles in [0,L]
    
    rho=0.*xg;
    for is=1:ns
        % RK-2, 1st step
        zps=squeeze(zp(is,:,:));
        pgmat=pginterp(zps,ng,dx,np,L);
        rhos=qs(is)*pgmat'*zps(:,4)/rhozero(is);
        rho=rho+rhos;
        if(is==1)
            pgmat1=pgmat;
        else
            pgmat2=pgmat;
        end
    end

    if(nfilter==1)
        rhof=fft(rho)/ng;
        af1=real(rhof(2)+rhof(ng)); bf1=imag(-rhof(2)+rhof(ng));
        if(only_ion==1) % 16-10-02
            % rhog=(af1*cos(2*pi/L*xg)+bf1*sin(2*pi/L*xg));
            % phig=rhog/(tau-k^2);
            % phig=rhog/(1+tau-Gamma0);
            Eg=-k*(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg))/(1+tau-Gamma0);
        else
            Eg=-(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg))/k;
        end
    elseif(nfilter==2) % 16-10-24 09:04, not ok
        rhof=fft(rho);
        ind=[2,ng]; % keep only the k=1*k mode
        rhof2=0.*rhof; rhof2(ind)=rhof(ind); % rho2=ng*ifft(rhof2); 
        Eg0=-ifft(rhof2./(1i*k));
        Eg=Eg0+conj(Eg0);% Eg=imag(Eg0);
    else
        Eg=cumsum(rho)*dx;
        Eg=Eg-mean(Eg);
    end

    for is=1:ns
        if(is==1)
           pgmat=pgmat1;
        else
           pgmat=pgmat2;
        end
        vdot1=pgmat*Eg;
        vdot2(1,:,:)=vdot1;
        vdot(1,:,:)=qs(is)/ms(is)*vdot1;
        wdot=(zp(is,:,3)-nonlinears(is)*zp(is,:,4)).*vdot2.*(kappan+...
            0.5*kappat*(-1.0+zp(is,:,2).^2/vts(is)^2)+zp(is,:,2)*theta*qs(is)/Ts(is));

        czp(is,:,1)=zp(is,:,1)+theta*zp(is,:,2)*dt/2;
        czp(is,:,2)=zp(is,:,2)+nonlinears(is)*theta*vdot*dt/2;
        czp(is,:,3)=zp(is,:,3);
        czp(is,:,4)=zp(is,:,4)+wdot*dt/2;
    end

    czp(:,:,1)=mod(czp(:,:,1)+10*L,L);
    rho=0.*xg;
    for is=1:ns
        % RK-2, 2nd step
        czps=squeeze(czp(is,:,:));
        pgmat=pginterp(czps,ng,dx,np,L);
        rhos=qs(is)*pgmat'*czps(:,4)/rhozero(is);
        rho=rho+rhos;
        if(is==1)
            pgmat1=pgmat;
        else
            pgmat2=pgmat;
        end
    end
    
    if(nfilter==1)
      %  rhof=fft(rhos)/ng; % 16-09-18 00:00 wrong! found by ZhangY
        rhof=fft(rho)/ng;
        af1=real(rhof(2)+rhof(ng)); bf1=imag(-rhof(2)+rhof(ng));
        if(only_ion==1) % 16-10-02
            Eg=-k*(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg))/(1+tau-Gamma0);
        else
            Eg=-(-af1*sin(2*pi/L*xg)+bf1*cos(2*pi/L*xg))/k;
        end
    elseif(nfilter==2) % 16-10-24 09:04
        rhof=fft(rho);
        ind=[2,ng]; % keep only the k=1*k mode
        rhof2=0.*rhof; rhof2(ind)=rhof(ind); % rho2=ng*ifft(rhof2);        
        Eg0=-ifft(rhof2./(1i*k));
        Eg=Eg0+conj(Eg0);% Eg=imag(Eg0);
    else
        Eg=cumsum(rho)*dx;
        Eg=Eg-mean(Eg);
    end
    
    Pt0(1:2)=0;
    for is=1:ns
        if(is==1)
           pgmat=pgmat1;
        else
           pgmat=pgmat2;
        end
        vdot1=pgmat*Eg;
        vdot2(1,:,:)=vdot1;
        vdot(1,:,:)=qs(is)/ms(is)*vdot1;
        wdot=(czp(is,:,3)-nonlinears(is)*czp(is,:,4)).*vdot2.*(kappan+...
            0.5*kappat*(-1.0+czp(is,:,2).^2/vts(is)^2)+czp(is,:,2)*theta*qs(is)/Ts(is));

        zp(is,:,1)=zp(is,:,1)+theta*czp(is,:,2)*dt;
        zp(is,:,2)=zp(is,:,2)+nonlinears(is)*theta*vdot*dt;
        zp(is,:,3)=zp(is,:,3);
        zp(is,:,4)=zp(is,:,4)+wdot*dt;
        
%         denflux(is,it)=mean(vdot1); % density flux, 16-10-03
%         engflux(is,it)=mean(mean((czp(is,:,2).^2/vts(is)^2).*vdot2)); % energy flux
        denflux(is,it)=mean(vdot2(1,:,:).*czp(is,:,4)); % density flux, 16-10-23 22:54
        engflux(is,it)=mean(mean((czp(is,:,2).^2/vts(is)^2).*vdot2(1,:,:).*czp(is,:,4))); % energy flux
        Pt(is,it)=ms(is)*sum(czp(is,:,2).*czp(is,:,4))/np;
        dPt(is,it)=-(Pt(is,it)-Pt0(is))/dt;
        Pt0(is)=Pt(is,it);
    end

%     EEk(it)=0.5*(sum(zp(1,:,2).^2))/np*L; % kinetic energy
%     EEk(it)=0.5*(sum(zp(1,:,2).^2)+sum(zp(2,:,2).^2))/np*L;
%     EEk(it)=0.5*(sum(zp(1,:,4).*zp(1,:,2).^2)*ms(1)+sum(zp(2,:,4).*zp(2,:,2).^2)*ms(2))/np*L;
    EEk(it)=0.5*(sum(zp(1,:,4).*zp(1,:,2).^2)*ms(1))/np;
    EEk2(it)=0.5*(sum(zp(2,:,4).*zp(2,:,2).^2)*ms(2))/np;
    EEk3(it)=0.5*(sum(zp(1,:,2).^2)*ms(1))/np;
    EEk4(it)=0.5*(sum(zp(2,:,2).^2)*ms(2))/np;
    EEf(it)=0.5*sum(Eg.^2)*dx; % potential energy
    if(nfilter~=1)
        Egf=fft(Eg)/ng;
        af1=real(Egf(2)+Egf(ng)); bf1=imag(-Egf(2)+Egf(ng));
    end
	Egt(it)=(af1+1i*bf1)/k;
    Egt0(it)=Eg(floor(ng/2));
    Dent(it)=sum(zp(2,:,4))/np;
    Dent1(it)=sum(zp(1,:,4))/np;
    
    % diagnosis
    zp(:,:,1)=mod(zp(:,:,1)+10*L,L);
    if(mod(it,100)==1)
    subplot(331);
    plot(xg,rho,'LineWidth',2);
    title(['\rho, it=',num2str(it),'/',num2str(ntb)]);
    xlim([min(xg),max(xg)]);
    subplot(332);
    plot(xg,Eg,'LineWidth',2);
    title('Eg');
    xlim([min(xg),max(xg)]);
    subplot(333);
%     plot(zp(1,:,2),zp(1,:,3),'.'); xlabel('v'); ylabel('p');
%     xlim([-1.1*vmax,1.1*vmax]);
    plot((1:it)*dt,real((Egt)),(1:it)*dt,imag((Egt)),':','Linewidth',2); 
    xlabel('t'); ylabel('real(Eg)');
    xlim([0,it*dt]);
    subplot(334);
    scatter(zp(1,:,1),zp(1,:,2),zp(1,:,3),zp(1,:,3)+zp(1,:,4));
    xlabel('x'); ylabel('v'); title('ion');
    ylim([-1.1*vmax,1.1*vmax]); xlim([min(xg),max(xg)]);
    subplot(335);
%     plot((1:it)*dt,real(log(EEf)),'Linewidth',2); 
    plot((1:it)*dt,real(log(abs(Egt))),'Linewidth',2); 
    xlabel('t'); ylabel('log(Egf)');
    xlim([0,it*dt]);
	if(only_ion~=1)
        subplot(336);
        scatter(zp(2,:,1),zp(2,:,2),zp(2,:,3),zp(2,:,4));
        ylim([-vts(2)*vmax,vts(2)*vmax]); xlim([min(xg),max(xg)]);
        xlabel('x'); ylabel('v'); title('electron');
        subplot(338);
        plot((1:it)*dt,real(denflux(1,:)),(1:it)*dt,real(denflux(2,:)),':'...
            ,(1:it)*dt,real(dPt(2,:)),'--','Linewidth',2); 
        xlabel('t'); ylabel('<\Gamma_s>');
        xlim([0,it*dt]);
        subplot(339);
        plot((1:it)*dt,real(engflux(1,:)),(1:it)*dt,real(engflux(2,:)),'Linewidth',2); 
        xlabel('t'); ylabel('<Q_s>');
        xlim([0,it*dt]);
	end
    subplot(337);
    plot((1:it)*dt,EEf,'b-',(1:it)*dt,-(EEk+EEk2),'g--',(1:it)*dt,Dent,'r:','LineWidth',2);
    xlabel('t'); legend('EEf','-EEk','density',2); legend('boxoff');
    xlim([0,it*dt]);
    
    drawnow;
    end
end
runtime=runtime+cputime-runtimea;
%%
str=['gkpic1d_df_k=',num2str(k),',ng=',num2str(ng),'np=',...
    num2str(np),',nt=',num2str(ntb),',dt=',num2str(dt),...
    ',nl=',num2str(nonlinears(1)),',nfilter=',num2str(nfilter),...
    ',nload=',num2str(nload),',noelectron=',num2str(only_ion)];

print(gcf,'-dpng',[str,'_snap.png']);
%%
% close all;
h = figure('Unit','Normalized','position',...
    [0.02 0.1 0.5 0.75],'DefaultAxesFontSize',15);
t=(1:ntb)*dt;
EEf1=EEf;
subplot(221);
plot(t,EEf1,'b-',t,1*(-EEk2-EEk),'g--',t,Dent,'r:','LineWidth',2);
% plot(t,EEf1,'b-',t,(-4*EEk),'g--',t,-2*Dent1,'r:','LineWidth',2);
% plot(t,EEf1,'b-',t,real(-dPt(1,:)/sqrt(ms(1))+dPt(2,:)/sqrt(ms(2))),'g--',t,Dent,'r:','LineWidth',2);
title(['(a) k=',num2str(k),', np=',num2str(np)]);
xlabel('t');axis tight;

subplot(222);
% Find the corresponding indexes of the extreme max values
% lndE=log(sqrt(real((EEf(1:nt))))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
lndE=log(real((Egt(1:ntb))));
lndEa=log(abs((Egt(1:ntb))));
it0=floor(ntb*2/20); it1=floor(ntb*6.9/20);
yy=real(lndE(it0:it1));
plot(t,lndE,t,lndEa,'LineWidth',2); hold on;
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
subplot(223);
% plot(t,real(denflux(2,:)),t,real(dPt(2,:)/vts(2)),'--','Linewidth',2);
plot(t,real(denflux(2,:)),t,real(dPt(2,:)/5),'--','Linewidth',2);
title(['(c) runtime=',num2str(runtime),'s']);
axis tight;
subplot(224);
plot(t,real(Pt(1,:)),t,real(Pt(2,:)),'--','Linewidth',2);
axis tight;

% str=['gkpic1d_df_k=',num2str(k),',ng=',num2str(ng),'np=',...
%     num2str(np),',nt=',num2str(ntb),',dt=',num2str(dt),...
%     ',nl=',num2str(nonlinear),',nfilter=',num2str(nfilter)];
print(gcf,'-dpng',[str,'_history.png']);
