% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2013-06-03 14:53
% Code for solve 1D shock wave tube, MacCormack.
% function shock1d_mac
close all; clear; clc;
% global gama nx Lx dx rho p u E U;
gama=1.4; 
nx=1000; Lx=2; 
x=linspace(0,Lx,nx+1); dx=x(2)-x(1);

% init
rhoa=1.0; rhob=0.125;
ua=0; ub=0;
pa=1.0; pb=0.1;

ind=find(x>Lx/2);
rho=rhoa.*ones(1,nx+1); rho(ind)=rhob;
u=ua.*ones(1,nx+1); u(ind)=ub;
p=pa.*ones(1,nx+1); p(ind)=pb;

U(1,:)=rho;
U(2,:)=rho.*u;
U(3,:)=p./(gama-1)+0.5.*rho.*u.*u;
Ef=0.*U; E=0.*U; Uf=0.*U;

nt=100; TT=0.4;
T=0; it=0;

figure('Unit','Normalized','Position',[0.01 0.2 0.5 0.5]);
set(gcf,'DefaultAxesFontSize',15);
% for it=1:nt
while(T<=TT)
    
    % CFL
    vel=sqrt(gama.*p./rho)+abs(u);
    CFL=real(dx/max(vel));
    dt=0.8*CFL;
    
    T=T+dt;
    it=it+1;
    
    r=dt/dx;
    dnu=0.35;
    
    % ¿ª¹Øº¯Êý
    q=abs(abs(U(1,3:nx+1)-U(1,2:nx))-abs(U(1,2:nx)-...
        U(1,1:nx-1)))./abs(abs(U(1,3:nx+1)-U(1,2:nx))+...
        abs(U(1,2:nx)-U(1,1:nx-1))+1e-10);
    q=repmat(q,3,1);
    % artificial viscous
    Ef(:,2:nx)=U(:,2:nx)+0.5*dnu*q.*(U(:,3:nx+1)-2*U(:,2:nx)+U(:,1:nx-1));
    
    U(:,2:nx)=Ef(:,2:nx);
    
    im=1; ip=nx+1;
%     Ef=U2E(U,im,ip);
    rho(im:ip)=U(1,im:ip);
    u(im:ip)=U(2,im:ip)./rho(im:ip);
    p(im:ip)=(gama-1).*(U(3,im:ip)-0.5*rho(im:ip).*u(im:ip).*u(im:ip));
    Ef(1,im:ip)=U(2,im:ip);
    Ef(2,im:ip)=rho(im:ip).*u(im:ip).*u(im:ip)+p(im:ip);
    Ef(3,im:ip)=(U(3,im:ip)+p(im:ip)).*u(im:ip);
    
    % U(n+1/2,i+1/2)
    Uf(:,1:nx)=U(:,1:nx)-r.*(Ef(:,2:nx+1)-Ef(:,1:nx));  
    
    im=1; ip=nx;  
    % E(n+1/2,i+1/2)
%     Ef=U2E(Uf,im,ip);
    rho(im:ip)=Uf(1,im:ip);
    u(im:ip)=Uf(2,im:ip)./rho(im:ip);
    p(im:ip)=(gama-1).*(Uf(3,im:ip)-0.5*rho(im:ip).*u(im:ip).*u(im:ip));
    Ef(1,im:ip)=Uf(2,im:ip);
    Ef(2,im:ip)=rho(im:ip).*u(im:ip).*u(im:ip)+p(im:ip);
    Ef(3,im:ip)=(Uf(3,im:ip)+p(im:ip)).*u(im:ip);
    
    % U(n+1,i)
    U(:,2:nx)=0.5*(U(:,2:nx)+Uf(:,2:nx))-0.5*r.*(Ef(:,2:nx)-Ef(:,1:nx-1));
    
    % b.c.
    U(:,1)=U(:,2);
    U(:,nx+1)=U(:,nx);
    
    % plot
%     if(it==nt || mod(it,floor(nt)/10)==1)
    if(mod(it,100)==1)
        plot(x,real(rho),'g',x,real(u),'r:',x,real(p),'Linewidth',2);
        legend('\rho','u','p');
        xlabel('x');xlim([min(x),max(x)]);ylim([-0.1,1.5]);
        title(['T=',num2str(T)]);
        pause(0.1);
    end

end
% end
% 
% function E=U2E(U,im,ip)
%     global gama rho p u;
%     rho(im:ip)=U(1,im:ip);
%     u(im:ip)=U(2,im:ip)./rho(im:ip);
%     p(im:ip)=(gama-1).*(U(3,im:ip)-0.5*rho(im:ip).*u(im:ip).*u(im:ip));
%     E(1,im:ip)=U(2,im:ip);
%     E(2,im:ip)=rho(im:ip).*u(im:ip).*u(im:ip)+p(im:ip);
%     E(3,im:ip)=(U(3,im:ip)+p(im:ip)).*u(im:ip);
% end


