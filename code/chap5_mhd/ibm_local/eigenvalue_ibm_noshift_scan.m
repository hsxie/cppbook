% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2011-10-30 18:25
% Ideal ballooning mode, eigenvalue solver, using matrix method
% Correction: 2012-02-18 20:01
% Update: using sparse() and eigs(), 2012-12-26 11:41
% 2013-03-07 09:47 without shift
close all;clear;clc;
Nx=1024*32;xmin=0;xmax=40;
% M=zeros(Nx);
dx=(xmax-xmin)/Nx;
dx2=dx*dx;
x=(xmin+dx):dx:xmax;

na=41;ns=41;
amin=0.0;amax=2.0;smin=0.0;smax=2.0;
aa=amin:(amax-amin)/(na-1):amax;
ss=smin:(smax-smin)/(ns-1):smax;
[AA,SS]=meshgrid(aa,ss);
WW=0.*AA;
runtime=cputime;
for ia=1:na
    for is=1:ns
        a=aa(ia);s=ss(is); % alpha, s
        
%         p=2.0.*(s.*x-a.*sin(x)).*(s-a.*cos(x));
%         q=a.*(cos(x)+(s.*x-a.*sin(x)).*sin(x));
%         r=-(1.0+(s.*x-a.*sin(x)).^2);
        
        p=2.0.*(s.*x).*(s);  % without shift, drop a.*sin(x) terms
        q=a.*(cos(x)+(s.*x).*sin(x));
        r=-(1.0+(s.*x).^2);

        
        I=zeros(3*Nx-2,1); J=I; S=I;
        I(1)=1; J(1)=1; S(1)=1.0/(dx2)-p(1)/(2.0*dx*r(1))+q(1)/r(1);
        for j=2:Nx
           I(j)=j; J(j)=j;
           S(j)=2.0/(dx2)+q(j)/r(j);
           I(Nx+j-1)=j; J(Nx+j-1)=j-1;
           S(Nx+j-1)=-1.0/(dx2)-p(j)/(2.0*dx*r(j));
           I(2*Nx+j-2)=j-1; J(2*Nx+j-2)=j;
           S(2*Nx+j-2)=-1.0/(dx2)+p(j-1)/(2.0*dx*r(j-1));
        end

        MS=sparse(I,J,S,Nx,Nx);

        sm=-2; % set closest parameter, default eigs 'sm' is 0, not enough
%         [V,D]=eigs(MS,2,sm);
        w2=eigs(MS,2,sm); % (omega/omega_A)^2
        w=sqrt(w2);
        
%         M(1,1)=1.0/(dx2)-p(1)/(2.0*dx*r(1))+q(1)/r(1);
%         for j=2:Nx
%             M(j,j)=2.0/(dx2)+q(j)/r(j);
%             M(j-1,j)=-1.0/(dx2)+p(j-1)/(2.0*dx*r(j-1));
%             M(j,j-1)=-1.0/(dx2)-p(j)/(2.0*dx*r(j));
%         end
% 
%         [V,D]=eig(M);
% 
%         w2=eig(M); % (omega/omega_A)^2
%         w=sqrt(w2);

        gammamax=max(imag(w));
        ind=find(imag(w)==gammamax);
        omega=w(ind(1));
        
        WW(is,ia)=omega; % not WW(ia,is)!!
    end
end
runtime=cputime-runtime;

set(gcf,'DefaultAxesFontSize',15);

subplot(221);mesh(AA,SS,imag(WW));xlabel('\alpha');ylabel('s');
zlabel('\gamma_{MHD}/\omega_A');title('Stable domain in s-\alpha space');
subplot(222); % contour(AA,SS,imag(WW));
pcolor(AA,SS,imag(WW)); shading interp;
xlabel('\alpha');ylabel('s');title(['Run time:',num2str(runtime),'s']);

row=floor(0.8/smax*(ns-1))+1;
subplot(223);plot(AA(row,:),imag(WW(row,:)),'LineWidth',2);xlabel('\alpha');
ylabel('\gamma_{MHD}/\omega_A');title(['s=',num2str(SS(row,1))]);

col=floor(0.8/amax*(na-1))+1;
subplot(224);plot(SS(:,col),imag(WW(:,col)),'LineWidth',2);xlabel('s');
ylabel('\gamma_{MHD}/\omega_A');title(['\alpha=',num2str(AA(1,col))]);

print(gcf,'-dpng',['stable_domain_noshift',',xmax=',num2str(xmax),...
    ',nx=',num2str(Nx),'.png']);
save s_a_noshift_scan;
