% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2016-04-29 22:39
% IVP to solve 1D kinetic ITG dispersion relation, d_t g(theta,v_par,v_perp)
% Benchmark with Dong1992PoF HD7
% Rewrite 2017-02-06 21:33
% 22:27 not ok
% 17-02-10 17:49 seems can agree with Dong92
% 20:39 optimize the bessel to speed up, seems faster 4x, e.g., 79s -> 22s
% typical case 128*32*32*5000 grids, cputime ~ 1200s
% 17-02-13 16:51 try RK2 -> RK4
% 23:44 fixed a bug on RK4 3rd step 0.5*dt should be dt
% 17-02-16 14:06 slightly update

close all;clear;clc;
irun=0; %
nt=2000;%
runtime=cputime;
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

    nvx=2*32/1; nvy=2*32/1; nx=128/1; dt=0.0025; % parameters
    xmax=4*pi; xmin=-xmax;
    vxmax=5.0; vxmin=-vxmax;
    vymax=5.0; vymin=0.0;

    dvx=(vxmax-vxmin)/nvx;
    dvy=(vymax-vymin)/nvy;
    dx=(xmax-xmin)/nx; % theta grid
    d2x=2*dx; d2x_inv=1.0/d2x;

    Vx=vxmin:dvx:vxmax; Vy=vymin:dvy:vymax; X=xmin:dx:xmax;  
%     Vy=Vy+dvy/2; % vy better start at dvy/2? 17-02-06 21:56
    [x,vx,vy]=ndgrid(X,Vx,Vy); % vx and vy grids

    F0=exp(-(vx.^2+vy.^2)/2.0)/sqrt((2*pi)^3); % initial distribution
    g=0.000.*F0;

    phi=0.1*exp(-X.^2/1.0^2);
    phit=zeros(nt,1); phit1=phit; phit2=phit;
    
    cg=g;
    cphi=phi;
    
    dg1=0.*g;dg2=0.*g;dg3=0.*g;dg4=0.*g;

    wdi=-kt;
    wsi=-kt/epsn;
    qR_inv=1/q;

    wD=wdi*(cos(x)+s*(x-tk).*sin(x)).*(vx.^2+vy.^2/2);
    wT=wsi*(1.0+(0.5*(vx.^2+vy.^2)-1.5)*etai);
    
    % 17-02-10 20:40, optimize, calculate them before time push
    ky=sqrt(kt^2.*(1+s^2*(x-tk).^2));
    ky0=sqrt(kt^2*(1+s^2*(X-tk).^2));
    Gamma0=besseli(0,ky0.^2,1);
    J0=besselj(0,ky.*vy);
    dJ0=-besselj(1,ky.*vy).*vy./ky.*(x-tk)*(kt*s)^2;
    J0FM=J0.*F0;
    dJ0FM=dJ0.*F0;
    hcoef=1./(1+1/tau-Gamma0);

    nta=1; ntb=nt;
else
    nta=ntb+1; ntb=nta+nt;
end

for it=nta:ntb
    
    % boundary in V space g must be null
%     g(nx+1,:,:)=g(nx,:,:); g(1,:,:)=g(2,:,:);
    g(nx+1,:,:)=0; g(1,:,:)=0; % 17-02-13 23:17
%     g(:,:,nvy+1)=0; g(:,nvx+1,:)=0;

    % RK-4, 1st step
    for ix=2:nx
        for ivx=1:nvx+1
            for ivy=1:nvy+1
                dg1(ix,ivx,ivy)=-vx(ix,ivx,ivy)*qR_inv*( ((g(ix+1,ivx,ivy)-g(ix-1,ivx,ivy))+...
                    J0FM(ix,ivx,ivy)*(phi(ix+1)-phi(ix-1)))*d2x_inv ...
                    +dJ0FM(ix,ivx,ivy)*phi(ix) )- 1i*wD(ix,ivx,ivy)*g(ix,ivx,ivy)+...
                    +1i*(wT(ix,ivx,ivy)-wD(ix,ivx,ivy))*J0FM(ix,ivx,ivy)*phi(ix);
            end
        end
    end
	cg=g+dg1*dt*0.5;
    g_int=(2*pi*sum(sum(cg.*vy.*J0,3),2)*dvx*dvy).';
    cphi=g_int.*hcoef;
    cphi(1)=0.0; cphi(nx+1)=0.0; cphi=smooth(cphi); % 16-05-24 20:41
    
    % RK-4, 2nd step
    for ix=2:nx
        for ivx=1:nvx+1
            for ivy=1:nvy+1
                dg2(ix,ivx,ivy)=-vx(ix,ivx,ivy)*qR_inv*( ((cg(ix+1,ivx,ivy)-cg(ix-1,ivx,ivy))+...
                    J0FM(ix,ivx,ivy)*(cphi(ix+1)-cphi(ix-1)))*d2x_inv ...
                    +dJ0FM(ix,ivx,ivy)*cphi(ix) )- 1i*wD(ix,ivx,ivy)*cg(ix,ivx,ivy)+...
                    +1i*(wT(ix,ivx,ivy)-wD(ix,ivx,ivy))*J0FM(ix,ivx,ivy)*cphi(ix);
                cg(ix,ivx,ivy)=g(ix,ivx,ivy)+dg2(ix,ivx,ivy)*dt*0.5;
            end
        end
    end
	cg=g+dg2*dt*0.5;
    g_int=(2*pi*sum(sum(cg.*vy.*J0,3),2)*dvx*dvy).';
    cphi=g_int.*hcoef;
    cphi(1)=0.0; cphi(nx+1)=0.0; cphi=smooth(cphi);
    
    % RK-4, 3rd step
    for ix=2:nx
        for ivx=1:nvx+1
            for ivy=1:nvy+1
                dg3(ix,ivx,ivy)=-vx(ix,ivx,ivy)*qR_inv*( ((cg(ix+1,ivx,ivy)-cg(ix-1,ivx,ivy))+...
                    J0FM(ix,ivx,ivy)*(cphi(ix+1)-cphi(ix-1)))*d2x_inv ...
                    +dJ0FM(ix,ivx,ivy)*cphi(ix) )- 1i*wD(ix,ivx,ivy)*cg(ix,ivx,ivy)+...
                    +1i*(wT(ix,ivx,ivy)-wD(ix,ivx,ivy))*J0FM(ix,ivx,ivy)*cphi(ix);
            end
        end
    end
% 	cg=g+dg3*dt*0.5; % wrong!!
	cg=g+dg3*dt; % 17-02-13 23:43
    g_int=(2*pi*sum(sum(cg.*vy.*J0,3),2)*dvx*dvy).';
    cphi=g_int.*hcoef;
    cphi(1)=0.0; cphi(nx+1)=0.0; cphi=smooth(cphi);
    
    % RK-4, 4th step
    for ix=2:nx
        for ivx=1:nvx+1
            for ivy=1:nvy+1
                dg4(ix,ivx,ivy)=-vx(ix,ivx,ivy)*qR_inv*( ((cg(ix+1,ivx,ivy)-cg(ix-1,ivx,ivy))+...
                    J0FM(ix,ivx,ivy)*(cphi(ix+1)-cphi(ix-1)))*d2x_inv ...
                    +dJ0FM(ix,ivx,ivy)*cphi(ix) )- 1i*wD(ix,ivx,ivy)*cg(ix,ivx,ivy)+...
                    +1i*(wT(ix,ivx,ivy)-wD(ix,ivx,ivy))*J0FM(ix,ivx,ivy)*cphi(ix);
            end
        end
    end
    g=g+(dg1+2*dg2+2*dg3+dg4)*dt/6.0;
            
    g_int=(2*pi*sum(sum(g.*vy.*J0,3),2)*dvx*dvy).';
    phi=g_int.*hcoef;
    phi(1)=0.0; phi(nx+1)=0.0; phi=smooth(phi);
        
	phit1(it)=mean(phi);
    phit2(it)=phi(floor(nx/2));
    if(it>=2 && mod(it,5)==0)
        subplot(221); plot(2:it,real(phit1(2:it)),2:it,...
            imag(phit1(2:it)),'r--','LineWidth',2);
        title(['phi, it=',num2str(it),'/',num2str(ntb)]);
        subplot(222); plot(X,real(phi),X,imag(phi),'Linewidth',2);
        subplot(223); plot(2:it,real(phit2(2:it)),2:it,...
            imag(phit2(2:it)),'r--','LineWidth',2);
        subplot(224); plot(2:it,log(real(phit1(2:it))),2:it,...
            log(imag(phit1(2:it))),'r--','LineWidth',2);
        drawnow;
    end
end
runtime=cputime-runtime;
phi0=phi;

%%
close all;
h = figure('Unit','Normalized','position',...
    [0.02 0.13 0.6 0.7],'DefaultAxesFontSize',15);

t=(2:ntb)*dt; phit=phit1(2:ntb);

subplot(221);
plot(t,real(phit),t,imag(phit),'r--','LineWidth',2); box on;
xlabel(['t, runtime=',num2str(runtime),'s']); legend('Re\phi','Im\phi',2); legend('boxoff');
title(['itg, k_\theta=',num2str(kt),', s=',num2str(s),', q=',num2str(q),...
    ', \epsilon_n=',num2str(epsn),', \eta_i=',num2str(etai),...
    ', \tau=',num2str(tau)]);

% Find the corresponding indexes of the extreme max values
% lndE=log(sqrt(real((phit(1:nt))))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
lndE=log(real((phit)));
it0=floor(ntb*4/20); it1=floor(ntb*19.9/20);
yy=lndE(it0:it1);
subplot(222);
plot(t,real(lndE),'LineWidth',2); box on;
hold on; axis tight;
extrMaxIndex = find(diff(sign(diff(yy)))==-2)+1;
t1=t(it0+extrMaxIndex(1));t2=t(it0+extrMaxIndex(end));
y1=yy(extrMaxIndex(1));y2=yy(extrMaxIndex(end));
plot([t1,t2],[y1,y2],'r*--','LineWidth',2);
nw=length(extrMaxIndex)-1;
omega=pi/((t2-t1)/nw);
gammas=(real(y2)-real(y1))/(t2-t1);
ws=omega+1i*gammas;
title(['(b) \omega^S=',num2str(ws,4),...
    ', nw=',num2str(nw)]);
xlabel(['\omega^T=',num2str(wt,3)]);

subplot(223);
% pcolor(vx(1,:,:),vy(1,:,:),real(g(1,:,:))); %log(abs(real(g)))
gx0=squeeze(g(floor(nx/2),:,:));
pcolor(squeeze(vx(1,:,:)),squeeze(vy(1,:,:)),real(gx0));
title(['vxmax=',num2str(vxmax),...
    ',vymax=',num2str(vymax),',nvx=',num2str(nvx),...
    ',nvy=',num2str(nvy),',nx=',num2str(nx)]);
subplot(224);
phi0=phi0./phi0(floor(nx/2));
plot(X,real(phi0),X,imag(phi0),'Linewidth',2);

w_ivp=ws;
th_ivp=X;
phi_ivp=phi0;
save('phi_th_ivp_c.mat','w_ivp','th_ivp','phi_ivp');

if(id==1)
run phi_hd7_dat; hold on;
plot(th_hd7,real(phi_hd7),'k:',th_hd7,imag(phi_hd7),'m:','Linewidth',2);
end

title(['xmax=',num2str(xmax),',dt=',num2str(dt),',ntb=',num2str(ntb)]);
print(gcf,'-dpng',['mgk1d_ivp_s=',num2str(s),',q=',num2str(q),...
    ',tau=',num2str(tau),',eta=',num2str(etai),',epsn=',num2str(epsn),...
    ',kt=',num2str(kt),',dt=',num2str(dt),',nx=',num2str(nx),...
    ',nvx=',num2str(nvx),',nvy=',num2str(nvy),...
    ',xmax=',num2str(xmax),',vxmax=',num2str(vxmax),...
    ',vymax=',num2str(vymax),'.png']);

