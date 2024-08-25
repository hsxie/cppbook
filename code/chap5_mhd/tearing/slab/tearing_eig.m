% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-06-22 11:03
% eigensolver for linear 2D slab tearing mode problem
% Ref: [1] FU Zhu-feng & HU You-qiu, 1995 book, p404-413
%      [2] Lee L.C. & Fu Z.F., 1986, JGR
% Test OK for growth rate, but seems not perfect for mode structure

close all; clear; clc;
runtime=cputime;

beta=0.2; alpha=0.3; Rm=8; gamma=3; % parameters
bc=1; % b.c.
eqm=0; % initial equilibrium type

l=1; L=10*l; dx=0.1*l;
x=-L:dx:L; d2x=2*dx; dx2=dx*dx;

% Index: -L --> j=0, L --> j=nx+1, matrix elements j=1, 2, ..., nx.
% Matlab array index 1:nx+2
nx=length(x)-2;

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


rhopj=(diff(rho(2:(nx+2)))+diff(rho(1:(nx+1))))./2/dx;
ppj=(diff(p(2:(nx+2)))+diff(p(1:(nx+1))))./2/dx;
Bzpj=(diff(Bz(2:(nx+2)))+diff(Bz(1:(nx+1))))./2/dx;

xj=x(2:(nx+1)); rhoj=rho(2:(nx+1)); pj=p(2:(nx+1)); Bzj=Bz(2:(nx+1));

% 1 -> rho1; 2 -> ux1; 3 -> uz1; 4 -> Bx1; 5 -> Bz1; 6 -> p1
% eigen matrix (6*nx)*(6*nx) dimensions
j=1:nx; jm=2:nx; jp=1:(nx-1);
% rho1
M12=sparse(j,j,-rhopj,nx,nx)+...
    sparse(jm,jp,rhoj(jm)/d2x,nx,nx)+...
    sparse(jp,jm,-rhoj(jp)./d2x,nx,nx);
M13=sparse(j,j,-1i*alpha*rhoj,nx,nx);

% ux1
M24=sparse(j,j,1i*alpha*Bzj./rhoj,nx,nx);
M25=sparse(j,j,-Bzpj./rhoj,nx,nx)+...
    sparse(jm,jp,Bzj(jm)./rhoj(jm)/d2x,nx,nx)+...
    sparse(jp,jm,-Bzj(jp)./rhoj(jp)/d2x,nx,nx);
M26=sparse(jm,jp,beta./(2*rhoj(jm))/d2x,nx,nx)+...
    sparse(jp,jm,-beta./(2*rhoj(jm))/d2x,nx,nx);

% uz1
M34=sparse(j,j,Bzpj./rhoj,nx,nx);
M36=sparse(j,j,-1i*alpha*beta./(2*rhoj),nx,nx);

% Bx1
M42=sparse(j,j,1i*alpha*Bzj,nx,nx);
M44=sparse(j,j,-(alpha^2+2/dx2)/Rm.*ones(1,nx),nx,nx)+...
    sparse(jm,jp,1/dx2/Rm.*ones(1,nx-1),nx,nx)+...
    sparse(jp,jm,1/dx2/Rm.*ones(1,nx-1),nx,nx);

% Bz1
M52=sparse(j,j,-Bzpj,nx,nx)+...
    sparse(jm,jp,Bzj(jm)/d2x,nx,nx)+...
    sparse(jp,jm,-Bzj(jp)/d2x,nx,nx);
M55=sparse(j,j,-(alpha^2+2/dx2)/Rm.*ones(1,nx),nx,nx)+...
    sparse(jm,jp,1/dx2/Rm.*ones(1,nx-1),nx,nx)+...
    sparse(jp,jm,1/dx2/Rm.*ones(1,nx-1),nx,nx);

% p1
M62=sparse(j,j,-ppj,nx,nx)+...
    sparse(jm,jp,gamma*pj(jm)/d2x,nx,nx)+...
    sparse(jp,jm,-gamma*pj(jp)/d2x,nx,nx);
M63=sparse(j,j,-1i*alpha*gamma*pj,nx,nx);

O=0.*M12;

% set matrix M
M=[O  M12 M13  O   O   O;
   O   O   O  M24 M25 M26;
   O   O   O  M34  O  M36;
   O  M42  O  M44  O   O ;
   O  M52  O   O  M55  O ;
   O  M62 M63  O   O   O];

% solve
sigma=0.1; % search eigenvalues aroud this value
[V,D]=eigs(M,6,sigma);
[row,col] = find(real(D)==max(max(real(D))));
w=1i*D(row,col);
rho1=V(j,col);
ux1=V(nx+j,col);
uz1=V(2*nx+j,col);
Bx1=V(3*nx+j,col);
Bz1=V(4*nx+j,col);
p1=V(5*nx+j,col);
    
h=figure('unit','normalized','Position',[0.01 0.1 0.65 0.7],...
    'DefaultAxesFontSize',15);
subplot(2,3,1);
plot(x,rho,'g',x,Bz,'b',x,p,'r--','LineWidth',2); 
xlabel('x'); xlim([-L,L]); axis tight;
legend('\rho_0','Bz_0','p_0',4); legend('boxoff');
title(['(a) \alpha=',num2str(alpha),', \beta=',...
    num2str(beta),', Rm=',num2str(Rm)]);

subplot(2,3,2);
if(nx<=2^9) % eig() supports only small dimensions, e.g., N<1000 
    FM=full(M);
    d=eig(FM); wtmp=1i*d;
    ind=find(imag(wtmp)>0);
    xmax=1.1*max(abs(real(wtmp))); 
    ymax=max(imag(wtmp)); ymin=min(imag(wtmp)); 
    plot(real(wtmp),imag(wtmp),'m.',real(wtmp(ind)),...
        imag(wtmp(ind)),'r+',[-xmax,xmax],[0,0],'g--','LineWidth',2);
    axis([-xmax,xmax 1.1*ymin-0.1*ymax 1.1*ymax-0.1*ymin]); 
end
xlabel('Re(\omega)'); ylabel('Im(\omega)');
title('(b) eigenvalues');

runtime=cputime-runtime;

subplot(2,3,3);
plot(xj,real(Bx1),xj,imag(Bx1),'r--','LineWidth',2); xlabel('x');
legend('Re(Bx_1)','Im(Bx_1)'); legend('boxoff');
title(['(c) nx=',num2str(nx),', runtime=',num2str(runtime,3),'s']);
subplot(2,3,4); 
plot(xj,real(ux1),xj,imag(ux1),'r--','LineWidth',2); xlabel('x');
legend('Re(ux1)','Im(ux1)'); legend('boxoff');
title(['(d) \omega=',num2str(w,3)]);
subplot(2,3,5);
plot(xj,real(p1),xj,imag(p1),'r--','LineWidth',2); xlabel('x');
legend('Re(p_1)','Im(p_1)'); legend('boxoff');
subplot(2,3,6);
plot(xj,real(Bz1),xj,imag(Bz1),'r--','LineWidth',2); xlabel('x');
legend('Re(Bz_1)','Im(Bz_1)'); legend('boxoff');

print('-dpng',['tearing_eig,alpha=',num2str(alpha),',beta=',...
    num2str(beta),',Rm=',num2str(Rm),',nx=',num2str(nx),'.png']);

