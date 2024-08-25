% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2016-04-06 15:36
close all; clear; clc;

M=1.0; c=1.0; e=1.0; m=1.0; mu=1.0;

% (r,theta,phi)->(x,y,z)
fx=@(r,theta,phi)r.*sin(theta).*cos(phi);
fy=@(r,theta,phi)r.*sin(theta).*sin(phi);
fz=@(r,theta,phi)r.*cos(theta);

% (chi,psi,zeta)->(r,theta,phi)
fr=@(chi,psi,zeta)M*sin(chi).^2./psi;
ftheta=@(chi,psi,zeta)chi;
fphi=@(chi,psi,zeta)zeta;

% B, B_chi, B_psi
fB=@(chi,psi,zeta)psi.^3.*sqrt(1+3*cos(chi).^2)./(M^2*sin(chi).^6);
fBchi=@(chi,psi,zeta)psi.^2.*(1+3*cos(chi).^2)./(M*sin(chi).^5);
fBpsi=@(chi,psi,zeta)-2*psi.*cos(chi)./(M*sin(chi).^4);

% d_psi B, d_chi B
fdchiB=@(chi,psi,zeta)-psi.^3*3.*cos(chi).*(3+...
    5*cos(chi).^2)./(M^2*sin(chi).^7.*sqrt(1+3*cos(chi).^2));
fdpsiB=@(chi,psi,zeta)3*psi.^2.*sqrt(1+3*cos(chi).^2)./(M^2*sin(chi).^6);


chi=pi/3; psi=5.0; zeta=0; vpar=0.1;
dt=0.0005; nt=1000; d=4; dt=dt/d; nt=nt*d;
zion=zeros(5,nt);

for it=1:nt
    
    zion(1,it)=chi; zion(2,it)=psi; zion(3,it)=zeta; zion(4,it)=vpar;
    zion(5,it)=dt*(it-1);
    
    % RK-2, 1st step
    rhs11=fB(chi,psi,zeta)/fBchi(chi,psi,zeta)*vpar;
    rhs21=0;
    rhs31=(c/e)*(mu+m*vpar^2/fB(chi,psi,zeta))*(fdpsiB(chi,psi,zeta)-...
        fBpsi(chi,psi,zeta)/fBchi(chi,psi,zeta)*fdchiB(chi,psi,zeta));
    rhs41=-mu/m*fB(chi,psi,zeta)/fBchi(chi,psi,zeta)*fdchiB(chi,psi,zeta);
    
    chi1=chi+dt*rhs11;
    psi1=psi+dt*rhs21;
    zeta1=zeta+dt*rhs31;
    vpar1=vpar+dt*rhs41;
    
    % RK-2, 2nd step
    rhs12=fB(chi1,psi1,zeta1)/fBchi(chi1,psi1,zeta1)*vpar1;
    rhs22=0;
    rhs32=(c/e)*(mu+m*vpar1^2/fB(chi1,psi1,zeta1))*(fdpsiB(chi1,psi1,zeta1)-...
        fBpsi(chi1,psi1,zeta1)/fBchi(chi1,psi1,zeta1)*fdchiB(chi1,psi1,zeta1));
    rhs42=-mu/m*fB(chi1,psi1,zeta1)/fBchi(chi1,psi1,zeta1)*fdchiB(chi1,psi1,zeta1);
    
    % RK-2, push
    chi=chi+0.5*dt*(rhs11+rhs12);
    psi=psi+0.5*dt*(rhs21+rhs22);
    zeta=zeta+0.5*dt*(rhs31+rhs32);
    vpar=vpar+0.5*dt*(rhs41+rhs42);
    
end
%%
r=fr(zion(1,:),zion(2,:),zion(3,:));
theta=ftheta(zion(1,:),zion(2,:),zion(3,:));
phi=fphi(zion(1,:),zion(2,:),zion(3,:));
x=fx(r,theta,phi);
y=fy(r,theta,phi);
z=fz(r,theta,phi);

zetatmp=zion(3,:);
ind=find(zetatmp<10*pi);
pt=ind(end); jt=1:pt;

close all;
h=figure('unit','normalized','Position',[0.01 0.17 0.67 0.55],...
    'DefaultAxesFontSize',15);
subplot(231);
plot(zion(5,jt),zion(1,jt),'-','Linewidth',2); 
xlabel('t'); ylabel('\chi'); axis tight;
subplot(232);
plot(zion(5,jt),zion(2,jt),'-','Linewidth',2); 
xlabel('t'); ylabel('\psi'); axis tight;
subplot(233);
plot(zion(5,:),zion(3,:),'Linewidth',2); 
xlabel('t'); ylabel('\zeta'); axis tight;
subplot(234);
plot(zion(5,jt),zion(4,jt),'Linewidth',2); 
xlabel('t'); ylabel('v_{||}'); axis tight;


subplot(235);
plot3(x,y,z,'Linewidth',2); box on;
xlabel('x'); ylabel('y'); zlabel('z');

subplot(236);
xlim([0,1]); ylim([0,1]); box off; axis off;
text(0,0.5,['chi=',num2str(zion(1,1)),...
    ', psi=',num2str(zion(2,1)),10,10,'zeta=',num2str(zion(3,1)),...
    ', vpar=',num2str(zion(4,1)),10,10,'dt=',num2str(dt),...
    ', nt=',num2str(nt)]);

print(gcf,'-dpng',['gc_orbit_chi=',num2str(zion(1,1)),...
    '_psi=',num2str(zion(2,1)),'_zeta=',num2str(zion(3,1)),...
    '_vpar=',num2str(zion(4,1)),'_dt=',num2str(dt),...
    '_nt=',num2str(nt),'.png']);
