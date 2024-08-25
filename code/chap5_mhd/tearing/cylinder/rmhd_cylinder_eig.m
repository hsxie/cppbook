% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-06-18 19:10
% Reduced MHD eigen solver, cylinder
% 2013-06-19 12:33 rewrite, seems OK

close all; clear; clc;

Nr=2^7;
r=linspace(0,1,Nr+2);
dr=r(2)-r(1);

% initial profiles
% q0=0.8; q1=-3.2; q2=4; 
% q0=1.1; q1=-1.8; q2=4;
q0=0.8; q1=0; q2=1; 
q=q0+q1*r+q2*r.^2;
qp=q1+2*q2.*r;
qpp=2*q2;
s=r.*qp./q;
sp=(r.*q.*qpp+q.*qp-qp.^2.*r)./q.^2;
eta=1e-6+0.*r;
nu=1e-6;
m=1; n=1;

rj=r(2:Nr+1); % define a temp radius parameter, r_j=j*dr, j=1,...,Nr
qj=q(2:Nr+1);
sj=s(2:Nr+1);
spj=sp(2:Nr+1);
etaj=eta(2:Nr+1);
rhojp=(1/dr^2+1./(2*rj*dr)); % rho_j^{+}
rhoj0=0.*rj-2/dr^2; % rho_j^{0}
rhojm=(1/dr^2-1./(2*rj*dr)); % rho_j^{-}


% set sparse matries DI, I, O, ...
j=1:Nr; jm=2:Nr; jp=1:(Nr-1);
DI = sparse(j,j,rhoj0,Nr,Nr)+sparse(jm,jp,rhojm(jm),Nr,Nr)+...
    sparse(jp,jm,rhojp(jp),Nr,Nr);

I=speye(Nr,Nr);
O=0.*DI;

tmp1=m^2./rj.^2;
Diag1 = sparse(j,j,tmp1,Nr,Nr);

tmp2=-m^2.*etaj./rj.^2;
Diag2 = sparse(j,j,tmp2,Nr,Nr);

tmp3=1i*(n-m./qj);
Diag3 = sparse(j,j,tmp3,Nr,Nr);

tmp4=-1i*((n-m./qj).*m^2./rj.^2+m./rj.*(spj./qj-sj.*(sj-2)./(rj.*qj)));
Diag4 = sparse(j,j,tmp4,Nr,Nr);

tmp5=m^4.*nu./rj.^4;
Diag5 = sparse(j,j,tmp5,Nr,Nr);

tmp6=etaj;
Diag6 = sparse(j,j,tmp6,Nr,Nr);

% set matrice A and B
B=[I O; O DI-Diag1];
C=[Diag2 Diag3; Diag4 Diag5];
D=[Diag6*DI O; Diag3*DI nu*(DI*DI-Diag1*DI-DI*Diag1)]; % '*' not '.*'!!
A=C+D;

% solve
sigma=0.1; % search eigenvalues aroud this value
[V,M]=eigs(A,B,6,sigma);
% [V,M]=eigs(A,B,6,'lr');
[row,col] = find(real(M)==max(max(real(M))));
w=1i*M(row,col);
psi=V(1:Nr,col);
phi=V((Nr+1):2*Nr,col);
norm=real(psi(find(abs(psi)==max(abs(psi)))));
psi=psi/norm; phi=phi/norm;
    
h=figure('unit','normalized','Position',[0.01 0.1 0.5 0.6],...
    'DefaultAxesFontSize',15);
subplot(2,2,1);
plot(r,q,r,s,[0,1],[1,1],'r--','LineWidth',2); 
xlabel('r'); xlim([0,1]); axis tight;
legend('q','s',2); legend('boxoff');
title(['(a) q=',num2str(q0),'+',num2str(q1),'r+',num2str(q2),...
    'r^2, \eta=',num2str(eta(1)),', \nu=',num2str(nu)]);

if(Nr<=2^9) % eig() supports only small dimensions, e.g., N<1000 
    FA=full(A); FB=full(B); FC=full(C); FD=full(D);
    d=eig(FA,FB); wtmp=1i*d;
    subplot(2,2,2);
    ind=find(imag(wtmp)>0);
    xmax=1.1*max(abs(real(wtmp))); 
    ymax=max(imag(wtmp)); ymin=min(imag(wtmp)); 
    plot(real(wtmp),imag(wtmp),'m.',real(wtmp(ind)),...
        imag(wtmp(ind)),'r+',[-xmax,xmax],[0,0],'g--','LineWidth',2);
    axis([-xmax,xmax 1.1*ymin-0.1*ymax 1.1*ymax-0.1*ymin]); 
    xlabel('Re(\omega)'); ylabel('Im(\omega)');
    title('(b) eigenvalues');
end

subplot(2,2,3);
plot(rj,real(psi),rj,imag(psi),'r--','LineWidth',2); xlabel('r');
legend('Re(\psi)','Im(\psi)'); legend('boxoff');
title(['(c) Nr=',num2str(Nr),', m=',num2str(m),', n=',num2str(n)]);
subplot(2,2,4); 
plot(rj,real(phi),rj,imag(phi),'r--','LineWidth',2); xlabel('r');
legend('Re(\phi)','Im(\phi)'); legend('boxoff');
title(['(d) \omega=',num2str(w,3)]);

print('-dpng',['Nr=',num2str(Nr),',m=',num2str(m),...
    ',n=',num2str(n),',q=',num2str(q0),'+',num2str(q1),...
    'r+',num2str(q2),'r^2.png']);

