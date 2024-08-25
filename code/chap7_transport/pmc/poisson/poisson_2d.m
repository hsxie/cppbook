% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2012-10-05 14:44
% rewrite from Rajeshwary Tayade, 2000
% To solve Laplace{u(x,y)} = g(x,y) using floating random walk
% define the function g(x,y) as 4, the 
% exact solution is u(x,y)=x^2+y^2
% boundary conditions: V(x,y) = x^2+y^2

close all; clear all; clc;
N=1e4; tol=1e-3;
Xlow=0; Xup=1;
Ylow=0; Yup=2;
[xx,yy]=meshgrid(Xlow:0.1*(Xup-Xlow):Xup,Ylow:0.1*(Yup-Ylow):Yup);
[nrow,ncol]=size(xx);
for irow=1:nrow
    for icol=1:ncol
        xstart=xx(irow,icol); ystart=yy(irow,icol);
        sum_u=zeros(1,N);
        Vb=zeros(1,N);
        uxy=0;
        % for each of the N walks
        for i=1:N
            x=xstart; % start at the same point
            y=ystart;
            rad=mindist(x,y,Xlow,Xup,Ylow,Yup);
            % each walk continues till u get a min distance close to a boundary
            while(rad>=tol)
                theta=rand*2*pi;
                x=x+rad*cos(theta);
                y=y+rad*sin(theta);
                % plot(x,y,'*-'); hold on; xlabel('x'); ylabel('y');
                % record g(x,y) *r*r
%                     sum_u(i)=sum_u(i)+(2*pi^2*sin(pi*x)*sin(pi*y)*rad*rad);
%                 sum_u(i)=sum_u(i)+(2*((1-6*x^2)*y^2*(1-y^2)+(1-6*y^2)*x^2*(1-x^2))*rad*rad);
                sum_u(i)=sum_u(i)+(-4*rad*rad);
                rad=mindist(x,y,Xlow,Xup,Ylow,Yup);
            end
            % determine which boundary is reached if Vb¡Ù0
            bnd= [x-Xlow, Xup-x, y-Ylow, Yup-y];
            [bd,j]=min(bnd);
    %         Vb(1,i)=(j==1)*(0)+(j==2)*(0)+(j==3)*(0)+(j==4)*(0);
%             Vb(1,i)=0;
            Vb(i)=x^2+y^2;
        end
        for i=1:N
            uxy=uxy+(Vb(i)+sum_u(i)/4);
        end
        uxy=uxy/N;
        uu(irow,icol)=uxy;
    end
end
figure; set(gcf,'DefaultAxesFontSize',15);
%  uue=sin(pi.*xx).*sin(pi.*yy);
% uue=-(xx.^2-xx.^4).*(yy.^2-yy.^4);
uue=xx.^2+yy.^2;
% minue=min(min(uue));maxue=max(max(uue));
subplot(121);surf(xx,yy,uu);
title('MC \phi{(x,y)}=x^2+y^2');xlabel('x');ylabel('y');
% zlim([1.1*minue-0.1*maxue,-0.1*minue+1.1*maxue]);
subplot(122);
% surf(xx,yy,uue);title('exact');xlabel('x');ylabel('y');
% zlim([1.1*minue-0.1*maxue,-0.2*minue+1.1*maxue]);
% figure;
midx=floor(nrow/2);
plot(xx(midx,:),uu(midx,:),'ro',xx(midx,:),uue(midx,:),':','LineWidth',2);
xlabel('x');ylabel(['\phi(x,',num2str(yy(midx,1)),')']);
title(['N=',num2str(N),' for per point']);
legend('MC','exact');legend('boxoff');
