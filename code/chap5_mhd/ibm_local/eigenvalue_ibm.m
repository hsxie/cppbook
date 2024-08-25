% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2011-10-30 18:25
% Ballooning mode, eigenvalue solver using matrix method
% Correction: 2012-02-18 13:17
% Update: 2012-12-26 09:56
% eig: very slow, and not accurate enough when N<256*8*n, n>2
% e.g., s=0.4, a=0.4, N=256*8, gamma=0.0793, run time > 1 min
% while shooting and initial give gamma ~ 0.091
% Rewrite from eig() to eigs() using sparse matrix
% N=1024*32*8, gamma=0.0911
close all;clear;clc;

s=0.4;a=0.8; % s, alpha
Nx=1024*32; xmin=0; xmax=40;
% M=zeros(Nx);
dx=(xmax-xmin)/Nx;
dx2=dx*dx;
x=(xmin+dx):dx:xmax;

p=2.0.*(s.*x-a.*sin(x)).*(s-a.*cos(x));
q=a.*(cos(x)+(s.*x-a.*sin(x)).*sin(x));
r=-(1.0+(s.*x-a.*sin(x)).^2);

% M(1,1)=1.0/(dx2)-p(1)/(2.0*dx*r(1))+q(1)/r(1);
% for j=2:Nx
%     M(j,j)=2.0/(dx2)+q(j)/r(j);
%     M(j-1,j)=-1.0/(dx2)+p(j-1)/(2.0*dx*r(j-1));
%     M(j,j-1)=-1.0/(dx2)-p(j)/(2.0*dx*r(j));
% end

I=zeros(3*Nx-2,1); J=I; S=I;
I(1)=1; J(1)=1; S(1)=1.0/(dx2)-p(1)/(2.0*dx*r(1))+q(1)/r(1);
for j=2:Nx
   I(j)=j; J(j)=j;S(j)=2.0/(dx2)+q(j)/r(j);
   I(Nx+j-1)=j;J(Nx+j-1)=j-1;S(Nx+j-1)=-1.0/(dx2)-p(j)/(2.0*dx*r(j));
   I(2*Nx+j-2)=j-1;J(2*Nx+j-2)=j;S(2*Nx+j-2)=-1.0/(dx2)+p(j-1)/(2.0*dx*r(j-1));
end

MS=sparse(I,J,S,Nx,Nx);
%%
set(gcf,'DefaultAxesFontSize',15);
sm=-5; % set closest parameter, default eigs 'sm' is 0, not enough
[V,D]=eigs(MS,2,sm);
w2=eigs(MS,20,sm); % (omega/omega_A)^2
w=sqrt(w2);
% hold on;
subplot(221);plot(real(w2),imag(w2),'x','LineWidth',2);xlabel('\lambda_r');
ylabel('\lambda_i');title('(a) eigenvalue of the matrix');

subplot(222);plot(real(w),imag(w),'x','LineWidth',2);xlabel('\omega');
ylabel('\gamma');title('(b) eigenfrequency');

% gammamax=max(imag(w));
% ind=find(imag(w)==gammamax);
eps=1e-3;
ind1=find(abs(imag(w2))<eps);
w2min=min(real(w2(ind1)));
ind=find(real(w2)==w2min);
omega=w(ind(1));
% Veff=a.*(cos(x)+(s.*x-a.*sin(x)).*sin(x))+omega.^2.*(1+(s.*x-a.*sin(x)).^2);

phi=V(1:Nx,ind(1))./V(1,ind(1));
subplot(212);
% plot(x,real(phi),'g',x,imag(phi),'g--',x,real(Veff),'r',x,imag(Veff),'r--
% ');
plot(x,real(phi),'b',x,imag(phi),'g--','LineWidth',2);
xlim([xmin,xmax]);ylim([-0.3,1.2]);
xlabel('\theta');ylabel('\phi');
title(['(c) max \gamma, ',...
    's= ',num2str(s),', \alpha= ',num2str(a),', \omega= ',num2str(omega)]);

print(gcf,'-dpng',['s=',num2str(s),',alpha=',num2str(a),',omega=',...
    num2str(omega),',xmax=',num2str(xmax),',nx=',num2str(Nx),'.png']);
% close all;

