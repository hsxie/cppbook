% Hua-sheng XIE, 2015-05-10
close all; clear; clc;
a=0.5; nt=2000; dt=0.02/2; dx=0.1; L=30.0; x=0.0:dx:L; nj=length(x)-1; 
figure('DefaultAxesFontSize',15);
method=2; 
for icase=1:2
    u=0.*x;
    if(icase==1)
        u(1:(nj+1))=exp(-(x-L/4).^2/1.8); subplot(211);        
    else
        u(floor(0.2*nj):floor(0.3*nj))=1.0; subplot(212);        
    end        
    x0=x; u0=u; j=2:nj;
    for it=1:nt
        if(method==1) % central difference 
            u(j)=u(j)-a*(u(j+1)-u(j-1))*dt/(2*dx);
        else % upwind
            u(j)=u(j)-a*(u(j)-u(j-1))*dt/(dx);
        end
       % u(1)=u(nj); u(nj+1)=u(2); % periodic b.c.
        if(mod(it,floor(nt/10))==1)
            plot(x0,u0,'k:',x,u,'b',x0+a*it*dt,u0,'r--','Linewidth',2); 
            legend('u_0(x)','u','u_0(x-at)'); legend('boxoff'); box on;
            xlim([min(x),max(x)]); ylim([-0.1,1.5]);  pause(0.1);
        end
    end
end