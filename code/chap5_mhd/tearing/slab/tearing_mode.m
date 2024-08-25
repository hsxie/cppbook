% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2012-11-21 20:19
% half spectral method for linear tearing mode problem
% Ref: [1] FU Zhu-feng & HU You-qiu, 1995 book, p404-413
%      [2] Lee L.C. & Fu Z.F., 1986, JGR
% 2013-01-12 15:03, fixed rhs and left b.c. bugs.
%    Test OK: both mode structure and linear growth rate.

close all; clear; clc;
runtime=cputime;

beta=0.2; alpha=0.3; Rm=8; gamma=3; % parameters
bc=1; % b.c.
eqm=0; % initial equilibrium type

l=1; dx=0.1*l; dt=0.001; nt=100000;
x=-10*l:dx:10*l; nx=length(x)-1;

% equilibrium
Bz0=1; RT=1; Bm=Bz0; pm=beta*Bm^2/2; rhom=pm/RT; p0=pm+Bm^2/2;
if(eqm==1)
    Bz=(abs(x)>=l).*Bz0.*x./abs(x+eps)+(abs(x)<l).*Bz0.*x/l;
    p=p0-Bz.^2/2; 
    rho=p./RT;
else
    Bz=Bz0*tanh(x./l);
    p=p0+sech(x./l).^2*Bz0^2/2; 
    rho=p./RT;
end

% drho=0.01*rand(1,nx+1);
% dux=0.01*rand(1,nx+1); 
% duz=0.01*rand(1,nx+1); 
% dBx=0.01*rand(1,nx+1); 
% dBz=0.01*rand(1,nx+1); 
% dp=0.01*rand(1,nx+1); 

drho=zeros(1,nx+1)+exp(-x.^2); 
% dux=0.0*drho; duz=-0.0*drho; dBx=0.0*drho+0.01;
% dBz=0.0*drho; dp=0.0*drho+0.01;
dux=0.01*drho; duz=-0.0*drho; dBx=0.01*drho; dBz=0.0*drho; dp=0.0*drho;

figure(2); set(gcf,'DefaultAxesFontSize',15);
for it=1:nt
    
   Am(it)=log(max(abs(dBx)));
   Amr(it)=log(max(abs(real(dBx))));
   t(it)=(it-1)*dt;
    
   drhotmp=drho;duxtmp=dux;duztmp=duz;
   dBxtmp=dBx;dBztmp=dBz;dptmp=dp;
   
   % r.h.s.
   drhodx=(diff(rho(2:nx+1))+diff(rho(1:nx)))./2/dx;
   dduxdx=(diff(dux(2:nx+1))+diff(dux(1:nx)))./2/dx;
   ddpdx=(diff(dp(2:nx+1))+diff(dp(1:nx)))./2/dx;
   dpdx=(diff(p(2:nx+1))+diff(p(1:nx)))./2/dx;
   ddBzdx=(diff(dBz(2:nx+1))+diff(dBz(1:nx)))./2/dx;
   dBzdx=(diff(Bz(2:nx+1))+diff(Bz(1:nx)))./2/dx;
   
   rhs1=-drhodx.*dux(2:nx)-rho(2:nx).*dduxdx...
       -1i*alpha.*rho(2:nx).*duz(2:nx); % version1 wrong: dBx --> duz !
   rhs2=-(beta/2)./rho(2:nx).*ddpdx-Bz(2:nx)./rho(2:nx).*ddBzdx-...
       dBzdx.*dBz(2:nx)./rho(2:nx)...
       +1i*alpha.*Bz(2:nx).*dBx(2:nx)./rho(2:nx);
   rhs3=-1i*alpha*(beta/2).*dp(2:nx)./rho(2:nx)+...
       dBzdx.*dBx(2:nx)./rho(2:nx); % version1 wrong: dBz --> dBx !
   
   d2dBxdx2=2*del2(dBx,dx); d2dBzdx2=2*del2(dBz,dx);
   
   rhs4=1i*alpha.*Bz(2:nx).*dux(2:nx)+d2dBxdx2(2:nx)/Rm-...
       alpha^2/Rm*dBx(2:nx);
   rhs5=-dBzdx.*dux(2:nx)-Bz(2:nx).*dduxdx+d2dBzdx2(2:nx)/Rm-...
       alpha^2/Rm*dBz(2:nx);
   rhs6=-dpdx.*dux(2:nx)-gamma*p(2:nx).*dduxdx-...
       1i*alpha*gamma.*p(2:nx).*duz(2:nx);
   
   drho(2:nx)=drho(2:nx)+rhs1*dt;
   dux(2:nx)=dux(2:nx)+rhs2*dt;
   duz(2:nx)=duz(2:nx)+rhs3*dt;
   dBx(2:nx)=dBx(2:nx)+rhs4*dt;
   dBz(2:nx)=dBz(2:nx)+rhs5*dt;
   dp(2:nx)=dp(2:nx)+rhs6*dt;
 
   if(bc==1) % b.c.: y'=0
       drho(1)=drho(2);
       dux(1)=dux(2);
       duz(1)=duz(2);
       dBx(1)=dBx(2);
       dBz(1)=dBz(2);
       dp(1)=dp(2);
           
       drho(nx+1)=drho(nx);
       dux(nx+1)=dux(nx);
       duz(nx+1)=duz(nx);
       dBx(nx+1)=dBx(nx);
       dBz(nx+1)=dBz(nx);
       dp(nx+1)=dp(nx);
   else % b.c.: y=0
       drho(1)=0;
       dux(1)=0;
       duz(1)=0;
       dBx(1)=0;
       dBz(1)=0;
       dp(1)=0;

       drho(nx+1)=0;
       dux(nx+1)=0;
       duz(nx+1)=0;
       dBx(nx+1)=0;
       dBz(nx+1)=0;
       dp(nx+1)=0;
   end
   
   if(it==nt || mod(it,floor(nt)/10)==1)
       subplot(221);plot(x,real(dBx),x,imag(dBx),'LineWidth',2);
       title('dBx-x');
       subplot(222);plot(x,real(dux),x,imag(dux),'LineWidth',2);
       title('dux-x');
       subplot(223);plot(x,real(dp),x,imag(dp),'LineWidth',2);
       title('dp-x');
       subplot(224);plot(t,Am,t,Amr,'LineWidth',2);
       title(['Am-t, t=',num2str((it-1)*dt)]);
       pause(0.05);
   end
end

runtime=cputime-runtime;
wi=(Am(nt)-Am(floor(nt/2)))/(t(nt)-t(floor(nt/2)));
subplot(224); axis tight; hold on;
title(['t=',num2str(nt*dt),', Am-t, \gamma=',num2str(wi)]);

print('-dpng',['tearing_mode_d,k=',num2str(alpha),',beta=',...
    num2str(beta),',Rm=',num2str(Rm),'.png']);

% plot eqm 
figure; set(gcf,'DefaultAxesFontSize',15);
subplot(221);plot(x,Bz,'LineWidth',2);title('Bz0-x');
subplot(222);plot(x,p,'LineWidth',2);title('p0-x');
subplot(223);plot(x,rho,'LineWidth',2);title('rho0-x');
subplot(224);title(['Run time=',num2str(runtime),'s']); axis off;
print('-dpng',['tearing_mode_eqm,k=',num2str(alpha),',beta=',...
    num2str(beta),',Rm=',num2str(Rm),'.png']);


