% Hua-shepsng XIE, huashepsngxie@gmail.com, FSC-ZJU, 2016-04-30 13:32
% matrix solve 1D kinetic ITG dispersion relation, g(theta,v_par,v_perp)
% Bepsnchmark with Dong1992PoF HD7
% 2016-09-03 23:52 Ref LiYY's version
% 16-09-06 21:38 benchmark with Dong92
% s = 1.0; kt = 0.45; tau = 1.0; etai= 2.5; epsn = 0.25; q = 1.0; gives 
% w_t=-0.607+0.258i;
% normalized by w=w*kt/w_{*e}
% Rewrite 2017-02-07 09:39
% 17-02-10 16:38 fixed a bug on RHS of omega*phi=..., miss a *hcoef
%   can almost agree with PIC now
% 17-02-16 15:06 seems still not well for w and g(x,vx,vy)

close all; clear; clc;
runtime=cputime;

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

vxmin = -5;
vxmax = -vxmin;
vymin = 0;
vymax = 5;
thmin = -4*pi;
thmax = -thmin;

nvx = 25; % change to hermite zeropoint?
nvy = 21; % change to laguerre zeropoint?
nth = 151;
totaln = nth*nvx*nvy;

dvx = (vxmax-vxmin)/(nvx-1);
dvy = (vymax-vymin)/(nvy-1);
dth = (thmax-thmin)/(nth-1);

vxc = linspace(vxmin,vxmax,nvx);
vyc = linspace(vymin,vymax,nvy);
thc = linspace(thmin,thmax,nth);

kt2=kt^2;
wdi=-kt;
wsi=-kt/epsn;
qR_inv=1/q;
fmc=1/sqrt((2*pi)^3);
ind = 1;

row = []; col = []; val = [];

for ith = 1:nth
    th = thc(ith);
    for ivx = 1:nvx
        vx = vxc(ivx);
        for ivy = 1:nvy
            vy = vyc(ivy);
            ky=sqrt(kt2*(1+s^2*(th-tk)^2));
            FM  = fmc*exp(-(vx^2+vy^2)/2); %
            J0 = besselj(0,ky*vy);
            dJ0=-besselj(1,ky*vy)*vy/ky*(th-tk)*(kt*s)^2;
            wd  = wdi*(cos(th)+s*(th-tk)*sin(th))*(vx^2+vy^2/2);
            ws  = wsi*(1+etai*(0.5*(vx^2+vy^2)-1.5));
        	rowind = (ith-1)*(nvx*nvy)+(ivx-1)*nvy+ivy;

            if (ith > 1)
            	row(ind) = rowind;
            	col(ind) = rowind - nvx*nvy;
            	val(ind) = 1i*qR_inv*vx/dth*2/3;

            	ind = ind + 1;

            	row(ind) = rowind;
            	col(ind) = totaln + ith-1;
            	val(ind) = 1i*qR_inv*vx*J0*FM/dth*2/3;

            	ind = ind + 1;
            end

            if (ith > 2)
                row(ind) = rowind;
            	col(ind) = rowind - nvx*nvy*2;
            	val(ind) = -1i*qR_inv*vx/dth/12;

            	ind = ind + 1;

            	row(ind) = rowind;
            	col(ind) = totaln + ith-2;
            	val(ind) = -1i*qR_inv*vx*J0*FM/dth/12;

            	ind = ind + 1;
            end

            row(ind) = rowind;
            col(ind) = rowind;
            val(ind) = wd;

            ind = ind + 1;

            row(ind) = rowind;
            col(ind) = totaln + ith;
            val(ind) = -(ws-wd)*J0*FM-1i*qR_inv*vx*dJ0*FM;

            ind = ind + 1;

            if (ith < nth)
            	row(ind) = rowind;
            	col(ind) = rowind + nvx*nvy;
            	val(ind) = -1i*qR_inv*vx/dth*2/3;

            	ind = ind + 1;

            	row(ind) = rowind;
            	col(ind) = totaln + ith+1;
            	val(ind) = -1i*qR_inv*vx*J0*FM/dth*2/3;

            	ind = ind + 1;
            end

            if (ith < nth-1)
            	row(ind) = rowind;
            	col(ind) = rowind + nvx*nvy*2;
            	val(ind) = 1i*qR_inv*vx/dth/12;

            	ind = ind + 1;

            	row(ind) = rowind;
            	col(ind) = totaln + ith+2;
            	val(ind) = 1i*qR_inv*vx*J0*FM/dth/12;

            	ind = ind + 1;
            end

        end
    end
end

for ith= 1:nth
    th = thc(ith);
    rowind = totaln + ith;
    bi=kt2*(1+s^2*(th-tk)^2);
    ky=sqrt(bi);
    fd=cos(th)+s*(th-tk)*sin(th);
    Gamma0 = besseli(0,bi,1);
    Gamma1 = besseli(1,bi,1);
    
    hcoef = 1/(1+1/tau-Gamma0);
    
    valphio=0;
    valphip=0;
    valphim=0;
    
    for ivx = 1:nvx
        vx = vxc(ivx);
        for ivy = 1:nvy
            vy = vyc(ivy);

            gind = (ith-1)*nvx*nvy+(ivx-1)*nvy+ivy;
            FM  = fmc*exp(-(vx^2+vy^2)/2);
            J0 = besselj(0,ky*vy);
            wd  = wdi*fd*(vx^2+vy^2/2);
            ws  = wsi*(1+etai*(0.5*(vx^2+vy^2)-1.5));

            if (ith > 1)
                row(ind) = rowind;
                col(ind) = gind-nvx*nvy;
                val(ind) = 2*pi*1i*qR_inv*vx*J0*vy/dth*2/3*dvx*dvy*hcoef;
                ind = ind+1;
            end
            
            if (ith > 2)
                row(ind) = rowind;
                col(ind) = gind-nvx*nvy*2;
                val(ind) = -2*pi*1i*qR_inv*vx*J0*vy/dth/12*dvx*dvy*hcoef;
                ind = ind+1;
            end

            row(ind) = rowind;
            col(ind) = gind;
            val(ind) = 2*pi*J0*wd*vy*dvx*dvy*hcoef;
            ind = ind+1;
                
            if (ith < nth)
                row(ind) = rowind;
                col(ind) = gind+nvx*nvy;
                val(ind) = -2*pi*1i*qR_inv*vx*J0*vy/dth*2/3*dvx*dvy*hcoef;
                ind = ind+1;
            end
            if (ith < nth-1)
                row(ind) = rowind;
                col(ind) = gind+nvx*nvy*2;
                val(ind) = 2*pi*1i*qR_inv*vx*J0*vy/dth/12*dvx*dvy*hcoef;
                ind = ind+1;
            end
            
        end
    end

    row(ind) = rowind;
    col(ind) = rowind;
%     val(ind) = wdi*fd*((2-bi)*Gamma0+bi*Gamma1)-wsi*((1-bi*etai)*Gamma0+bi*etai*Gamma1);
    val(ind) = (wdi*fd*((2-bi)*Gamma0+bi*Gamma1)-... % 17-02-10 16:36
        wsi*((1-bi*etai)*Gamma0+bi*etai*Gamma1))*hcoef;
    ind = ind + 1;
end

wholemat = sparse(row,col,val);

%%
M=wholemat;
wg=wt;
% wg=(-0.597+0.222i)*kt;
d=eigs(M,3,wg); %/kt
% d=eigs(M,2,'li');
[di,ind]=sort(imag(d),'descend');
d=d(ind);
[V,D]=eigs(M,1,d(1));

runtime=cputime-runtime;
%%
h = figure('Unit','Normalized','position',...
    [0.02 0.13 0.6 0.7],'DefaultAxesFontSize',15);

subplot(221);
plot(real(d),imag(d),'o','linewidth',2); hold on;
xlabel(['\omega_r, runtime=',num2str(runtime),'s']); ylabel('\omega_i');
title(['kt=',num2str(kt),',epsn=',num2str(epsn),',eta=',num2str(etai),...
    ',tau=',num2str(tau)]);

nd=nth*nvx*nvy; Nd=nd+nth;
jd=1;
V1=squeeze(V(:,jd));
phi=V1((nd+1):Nd);
% gxv=reshape(V1(1:nd),nth,nvx,nvy);
% gx0=squeeze(gxv(floor(nth/2),:,:));
gxv=reshape(V1(1:nd),nvy,nvx,nth); % 17-02-17 15:38
gx0=squeeze(gxv(:,:,floor(nth/2)));
gx0=gx0.';

plot(real(D(jd,jd)),imag(D(jd,jd)),'r*');

X=thc;
subplot(222);
phi=phi/phi(floor(nth/2));
plot(X,real(phi),X,imag(phi),'linewidth',2);
title(['\omega=',num2str(D(jd,jd),3)]);
xlabel(['\vartheta, \omega=',num2str(wt,3)]); ylabel('\phi');
if(id==1)
run phi_hd7_dat; hold on;
    plot(th_hd7,real(phi_hd7),'k:',...
        th_hd7,imag(phi_hd7),'m:','Linewidth',2);
end

th_matrix=thc;
phi_matrix=phi;

subplot(223);
[vx,vy]=meshgrid(vxc,vyc);
pcolor(vx',vy',real(gx0));
xlabel('v_{||}'); ylabel('v_\perp');
title(['vxmax=',num2str(vxmax),...
    ',vymax=',num2str(vymax),',nvx=',num2str(nvx),...
    ',nvy=',num2str(nvy),',nth=',num2str(nth)]);

subplot(224);
pcolor(vx',vy',imag(gx0));
xlabel('v_{||}'); ylabel('v_\perp');
title(['s=',num2str(s),',q=',num2str(q)]);

print(gcf,'-dpng',['mgk1d_matrix_s=',num2str(s),',q=',num2str(q),...
    ',tau=',num2str(tau),',eta=',num2str(etai),',epsn=',num2str(epsn),...
    ',kt=',num2str(kt),',nth=',num2str(nth),...
    ',nvx=',num2str(nvx),',nvy=',num2str(nvy),...
    ',thmax=',num2str(thmax),',vxmax=',num2str(vxmax),...
    ',vymax=',num2str(vymax),',w=',num2str(D(jd,jd)),'.png']);
