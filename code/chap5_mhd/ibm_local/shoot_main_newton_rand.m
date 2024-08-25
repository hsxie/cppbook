% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2012-02-19 10:29
% Ideal ballooning mode, shooting solver, using fun_newton_rand.m

close all;clear all;clc;

s=0.4;a=0.8; % s, alpha

nx=1001;xmin=0;xmax=40;dx=(xmax-xmin)/(nx-1);
tol=1e-3;
y0=[1,0];yn=0;

itmax=50;

found=0;
for II=1:1000
%     g(1)=0.0+0.9i;g(2)=0.0+0.65i; % initial guess (Q: How to set it for more cases?)

    g(1)=2*(rand(1)-0.5)+2*(rand(1)-0.5)*1i; % rand initial guess 
    g(2)=2*(rand(1)-0.5)+2*(rand(1)-0.5)*1i;
    
    for it=1:itmax
        w=g(it);
        [x,y]=ode45(@(x,y)fun_newton_rand(x,y,w,s,a),xmin:dx:xmax,y0);
        c(it)=y(end,1);
        if((abs(yn-c(it))<tol)&&(abs(y(end,2)-0)<tol)) break, end;
        if(it>=2)   % secant method
            g(it+1)=g(it)-(c(it)-yn)*(g(it)-g(it-1))/(c(it)-c(it-1));
        end
    end
    
    if(it<itmax)
        found=1;
        break
    end
end

set(gcf,'DefaultAxesFontSize',15);
plot(x,real(y(:,1)),x,imag(y(:,1)),'LineWidth',2);
xlim([xmin,xmax]);
if(found)
    title(['s=',num2str(s),',alpha=',num2str(a),',\omega=',num2str(w)],'fontsize',15);
    print(gcf,'-dpng',['s=',num2str(s),',alpha=',num2str(a),',omega=',num2str(w),'.png']);
else
    title('Havenot found a solution');
end
