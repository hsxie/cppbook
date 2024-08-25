% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-10-05 10:35
% Reduced MHD linear initial solver, cylinder
% 13:26 test OK, gamma_initial=0.0073414, gamma_eig=0.00734

close all; clear; clc;

nt=100000; dt=0.01;
Nr0=2^8; r=linspace(0,1,Nr0+2); dr=r(2)-r(1);

% initial profiles
% q0=0.8; q1=-3.2; q2=4;
% q0=1.1; q1=-1.8; q2=4;
% q0=0.8; q1=0; q2=1; 
q0=0.8; q1=0; q2=0.8; % J. McClenaghan GTC kink case
q=q0+q1*r+q2*r.^2;
qp=q1+2*q2.*r;
qpp=2*q2;
s=r.*qp./q;
sp=(r.*q.*qpp+q.*qp-qp.^2.*r)./q.^2;
eta=1e-6+0.*r;
nu=1e-6;
m=1; n=1;

Na=2; Nb=Nr0;
Nr=Nb+1-Na;

rj=r(Na:Nb); % define a temp radius parameter, r_j=j*dr, j=1,...,Nr
qj=q(Na:Nb);
sj=s(Na:Nb);
spj=sp(Na:Nb);
etaj=eta(Na:Nb);

rhojp=(1/dr^2+1./(2*rj*dr)); % rho_j^{+}
rhoj0=0.*rj-2/dr^2; % rho_j^{0}
rhojm=(1/dr^2-1./(2*rj*dr)); % rho_j^{-}

psij=0.00*exp(-(rj-0.2).^2/0.01);
uj=(1+0i)*0.01*exp(-(rj-0.3).^2/0.01);
phij=0.*rj;

% set sparse matries DI
j=1:Nr; jm=2:Nr; jp=1:(Nr-1);
DI = sparse(j,j,rhoj0-m^2./rj.^2,Nr,Nr)+sparse(jm,jp,rhojm(jm),Nr,Nr)+...
    sparse(jp,jm,rhojp(jp),Nr,Nr);

  
h=figure('unit','normalized','Position',[0.01 0.17 0.75 0.75],...
    'DefaultAxesFontSize',15);
subplot(331);
plot(r,q,r,s,[0,1],[1,1],'r--','LineWidth',2); 
xlabel('r'); xlim([0,1]); axis tight;
legend('q','s',2); legend('boxoff');
title(['q=',num2str(q0),'+',num2str(q1),'r+',num2str(q2),...
    'r^2, \eta=',num2str(eta(1)),', \nu=',num2str(nu)]);

gam(1)=1;
for it=1:nt
    
    Am(it)=log(max(abs(psij)));
    Amr(it)=log(max(abs(real(psij))));
    if(it>1)
        gam(it)=(Am(it)-Am(it-1))/dt;
    end
    t(it)=(it-1)*dt;

    phij=(DI\uj')';
    grad2perppsi=(DI*psij')';
    grad2perpu=(DI*uj')';
    
    rhs1=1i*(n-m./qj).*phij+etaj.*grad2perppsi;
    rhs2=1i*(n-m./qj).*grad2perppsi+nu*grad2perpu-...
        1i*m./rj.*(spj./qj-sj.*(sj-2)./(rj.*qj)).*psij;
    
    psij=psij+rhs1*dt;
    uj=uj+rhs2*dt;
    
   
   if(it==nt || mod(it,floor(nt)/10)==1)
       subplot(332);plot(rj,real(psij),rj,imag(psij),'r--','LineWidth',2);
       title('\Psi');
       subplot(333);plot(rj,real(phij),rj,imag(phij),'r--','LineWidth',2);
       title('\phi');
       subplot(334);plot(rj,real(uj),rj,imag(uj),'r--','LineWidth',2);
       title('U');
       subplot(335);plot(t,Am,t,Amr,'r--','LineWidth',2);
       title(['Am, t=',num2str((it-1)*dt)]); xlabel('t');
       subplot(336);plot(t,gam,'r',[0,t(it)],[0,0],'g--','LineWidth',2);
       title('\gamma'); xlabel('t');
       ylim([-1.5*abs(gam(it)),1.5*abs(gam(it))]);

       pause(0.02);
   end
   
end
%%
nta=floor(2*nt/3);
wi=(Am(nt)-Am(nta))/(t(nt)-t(nta));
subplot(335); axis tight; hold on;
title(['\Psi_{max}, t=',num2str(nt*dt),', \gamma=',num2str(wi)]);

%%
% map to 2D
the=0:pi/50:2*pi;
[rr,tt]=meshgrid(rj,the); % (r,theta)
[col,row]=size(rr);
psi=repmat(psij,col,1).*exp(-1i*m*tt);
phi=repmat(phij,col,1).*exp(-1i*m*tt);
U=repmat(uj,col,1).*exp(-1i*m*tt);
X=rr.*cos(tt); Z=rr.*sin(tt);
subplot(337); pcolor(X,Z,real(psi));
shading interp; axis equal; title('\psi, 2D contour');
colormap(hsv); colorbar('location','EastOutside');
subplot(338); pcolor(X,Z,real(phi));
shading interp; axis equal; title('\phi, 2D contour');
colormap(hsv); colorbar('location','EastOutside');
subplot(339); pcolor(X,Z,real(U));
shading interp; axis equal; title('{U}, 2D contour');
colormap(hsv); colorbar('location','EastOutside');

print('-dpng',['tearing_initial_Nr=',num2str(Nr),',m=',num2str(m),...
    ',n=',num2str(n),',q=',num2str(q0),'+',num2str(q1),...
    'r+',num2str(q2),'r^2.png']);


