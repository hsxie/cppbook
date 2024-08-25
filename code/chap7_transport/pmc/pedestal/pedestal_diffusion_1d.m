% heat.m
%
% Solution of the heat equation u_t=Du_xx
% using expectation of the heat kernel
% and random walk methods
%
close all;
clear all;
M=100;
D=1;
x=linspace(-5,5,M);
dx=x(2)-x(1);
ue=0.5*((x>-2.5)&(x<-0.5))+((x>0.5)&(x<2.5));
ue=ue/(dx*sum(ue));
N=12000;
xi=-2.5+2*rand(1,N/3);
xi(N/3+1:N)=0.5+2*rand(1,N-N/3);
xii=xi;
N2=round(N/M);
xi0=-2.5+2*rand(1,N2/3);
xi0(N2/3+1:N2)=0.5+2*rand(1,N2-N2/3);
m=1/(N*dx);
u=m*hist(xi,x);
plot(x,u,':.',x,ue);
dt=0.004;
lambda=D*dt/(dx)^2;
for t=1:100
    ue(2:M-1)=ue(2:M-1)*(1-2*lambda)+lambda*(ue(1:M-2)+ue(3:M));
    xi=xii+sqrt(2*D*dt*t)*randn(1,N);
    u=m*hist(xi,x);
    for i=1:M
        u2(i)=sum(exp(-(x(i)-xi0).^2./(4*D*dt*t)))*M/(N*sqrt(4*D*pi*dt*t));
    end
    plot(x,u,':.',x,ue,'-k',x,u2,':.r');
    axis([-5 5 0 0.6]);
    drawnow;
    e1(t)=sum((ue-u).^2);
    e2(t)=sum((ue-u2).^2);
end

figure;
semilogy((1:100)*dt,e1,(1:100)*dt,e2);

