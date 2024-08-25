% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2012-06-15 12:28
% Solving the Ex14.4 in Freidberg2007's book.
% Q: How to treat the boundary? A: dx should be large!
close all;clear;clc;
D=1.0; nt=1001; dt=0.0002; a=1.0; nj=50; dr=a/nj; r=0:dr:a;
% r=dr/2:dr:a+dr/2;
n=zeros(nt+1,nj+1);
n0=1.0;
n(1,:)=n0;
% n(1,:)=n0.*(1.0+0.99.*cos(pi.*r));
% n(1,:)=n0.*exp(-r.^2./(0.3*a)^2);
% n(1,nj+1)=0;
set(gcf,'DefaultAxesFontSize',15);
plot(r,n(1,:),'r:','LineWidth',2);hold on;
for it=1:nt
    if(mod(it,floor(nt/10))==1)        
        if(it<=1)
            plot(r,n(1,:),'r:','LineWidth',2);hold on;
        else
            plot(r,n(it,:),'g','LineWidth',2);
        end
        xlabel('r');ylabel('n');
        title(['1D cylinder diffusion, D=',num2str(D),', dt=',num2str(dt),...
            ', dr=',num2str(dr,3),', a=',num2str(a)]); ylim([0,1.5]);
        [ym,idx]=max(n(it,:),[],2);
        text(r(idx),ym,['t=',num2str((it-1)*dt)]);
    end
    
    for j=2:nj
        n(it+1,j)=n(it,j)+D*((n(it,j+1)+n(it,j-1)-2*n(it,j))/(dr^2)+(n(it,j+1)-n(it,j-1))/(2*r(j)*dr))*dt;
    end
    n(it+1,1)=n(it+1,2);
%     n(it+1,nj+1)=n(it+1,nj);
    n(it+1,nj+1)=0;
    pause(0.01);
end
print(gcf,'-dpng','diffusion_1d_r.png');
