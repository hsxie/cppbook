% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2012-06-15 12:07
% Solving the diffusion equation, dn/dt=D*d^2n/dx^2.
% Q: How to treat numerical diffusion? A: dx should be large!
close all;clear;clc;
D=1.0; nt=1001; dt=0.0001; L=1.0; nj=30; dx=L/nj; x=0:dx:L;
n=zeros(nt+1,nj+1);
n(1,:)=2.0.*sin(pi.*x./L)+1.0.*sin(3*pi.*x./L)+0.5.*sin(5*pi.*x./L);
set(gcf,'DefaultAxesFontSize',15);
for it=1:nt
    if(mod(it,floor(nt/10))==1)
        if(it<=1)
            plot(x,n(1,:),'r:','LineWidth',2);hold on;
        else
            plot(x,n(it,:),'g','LineWidth',2);
        end
        xlabel('x');ylabel('n');
        title(['1D diffusion, D=',num2str(D),', dt=',num2str(dt),...
            ', dx=',num2str(dx,3),', L=',num2str(L)]);%ylim([0,2.5]);
        [ym,idx]=max(n(it,:),[],2);
        text(x(idx),ym,['t=',num2str((it-1)*dt)]);
    end
    for j=2:nj
        n(it+1,j)=n(it,j)+D*(n(it,j+1)+n(it,j-1)-2*n(it,j))/(dx^2)*dt;
    end
    n(it+1,1)=0; n(it+1,nj+1)=0; pause(0.01);
end
print(gcf,'-dpng','diffusion_1d_x.png');
