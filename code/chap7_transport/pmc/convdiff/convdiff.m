% convdiff.m, 2012-10-03 20:49, huashengxie@gmail.com
% brown from Lorenzo Pareschi, rewritten by Hua-sheng XIE, IFTS-ZJU
%
% Solution of the convection-diffusion equation u_t+au_x=Du_xx
% using free particle transport and free boundary condition
%
close all; clear all; clc;

figure; set(gcf,'DefaultAxesFontSize',15);

M=200; Nt=100; a=1; D=0.1; N=1e5; dt=0.01;
x=linspace(-5,5,M);
dx=x(2)-x(1);
ue=0.5*((x>-4)&(x<-2))+((x>-0.5)&(x<0.5));
ue=ue/(dx*sum(ue));
xi=-4+2*rand(1,N/2);
xi(N/2+1:N)=-0.5+rand(1,N-N/2);

m=1/(N*dx);
u=m*hist(xi,x);

lambda2=D*dt/(dx)^2;
lambda1=a*dt/(2*dx);

subplot(121); plot(x,u,'r:o',x,ue,'-','LineWidth',2);
title(['a=',num2str(a),', D=',num2str(D),', dt=',num2str(dt)]);
text(-4,0.4,'t=0','FontSize',15);
axis([-5 5 0 0.6]); xlabel('x');ylabel('f');

for t=1:Nt
    ue(2:M-1)= lambda1*(ue(1:M-2)-ue(3:M))+(1-2*lambda2)*ue(2:M-1)+lambda2*(ue(1:M-2)+ue(3:M));
    xi=xi+a*dt+sqrt(2*D*dt)*randn(1,N);
    u=m*hist(xi,x);
    subplot(122);plot(x,u,'r:o',x,ue,'-','LineWidth',2);
    title(['N=',num2str(N),', M=',num2str(M)]);
    text(-4,0.4,['t=',num2str(dt*t)],'FontSize',15);
    axis([-5 5 0 0.6]);xlabel('x');ylabel('f');
    legend('MC','exact'); legend('boxoff');
    drawnow;
    pause(.01);
%     e1(t)=sum((ue-u).^2);
end
% figure;
% semilogy((1:Nt)*dt,e1);
