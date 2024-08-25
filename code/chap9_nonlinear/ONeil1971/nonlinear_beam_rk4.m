close all; clear;clc;
np=1000; dt=0.01; nt=1800; L=2.0*pi;
xp=linspace(0,L-L/np,np)'; % xpj(t=0)
phi=0.01;
phi_t=[];Ek_t=[];Omega_t=[];Gamma_t=[];

omega(1)=0.5*(-1+sqrt(3)*1i);
% position perurbation use the 1st order approximation (Zonca lecture 3:3-7)
xp=xp+2*real((0+1i)*phi(1)*exp((0+1i)*xp)/omega(1)^2);
vp=0.*xp+2*real(phi(1)*exp((0+1i)*xp)/omega(1));
sum(exp(-1i*xp))

% % bc
t=linspace(0,nt*dt,nt);
tmx=max(t);
ta=0; %ta=tmx*0.25;
h=figure('unit','normalized','position',[0.1,0.1,0.6,0.7],...
    'DefaultAxesFontSize',12);
jp=1;
for it=1:nt

    % RK-4, 1st step
    u1=vp;
    a1=2.0*real(-(1i*phi).*exp(1i.*xp));
    g1=-1i*sum(exp(-1i.*xp))/np;
    % RK-4, 2nd step
    xtmp=xp+0.5.*dt.*u1;
    vtmp=vp+0.5.*dt.*a1;
    phitmp=phi+0.5*dt*g1;
    u2=vtmp;
    a2=2.0*real(-(1i*phitmp).*exp(1i.*xtmp));
    g2=-1i*sum(exp(-1i.*xtmp))/np;
    % RK-4, 3rd step
    xtmp=xp+0.5.*dt.*u2;
    vtmp=vp+0.5.*dt.*a2;
    phitmp=phi+0.5*dt*g2;
    u3=vtmp;
    a3=2.0*real(-(1i*phitmp).*exp(1i.*xtmp));
    g3=-1i*sum(exp(-1i.*xtmp))/np;
    % RK-4, 4th step
    xtmp=xp+dt.*u3;
    vtmp=vp+dt.*a3;
    phitmp=phi+dt*g3;
    u4=vtmp;
    a4=2.0*real(-(1i*phitmp).*exp(1i.*xtmp));
    g4=-1i*sum(exp(-1i.*xtmp))/np;
    % RK-4, push
    xp=xp+dt./6.0.*(u1+2.0.*u2+2.0.*u3+u4);
    vp=vp+dt./6.0.*(a1+2.0.*a2+2.0.*a3+a4);
    phitmp=phi;
    phi=phi+dt/6.0*(g1+2.0*g2+2.0*g3+g4);
    % Omega
    Omega=1i*(phi-phitmp)/dt/phi;
    Gamma=(abs(phi)-abs(phitmp))/dt/abs(phi);
    
    % bc
    xp=xp./L+100.0;
    xp=L.*(xp-floor(xp));
    
    % plot
    if(mod(it,floor(nt/6))==1 && jp<=6)
        subplot(2,3,jp); jp=jp+1;        
        plot(xp,vp,'.');xlim([0,L]); ylim([-4.0,3.0]);
        text(0.1*L,2.0,num2str(it*dt,'t=%04.2f'));
        xlabel('POSITION(\xi)');ylabel('VELOCITY(\xi'')');
    end
    phi_t=[phi_t,phi]; % record phi(t)
    Ek_t=[Ek_t,mean(vp)];
    Omega_t=[Omega_t,Omega];
    Gamma_t=[Gamma_t,Gamma];
end
%%
h=figure('unit','normalized','position',[0.02,0.1,0.6,0.7],...
    'DefaultAxesFontSize',12);
subplot(221);plot(t,abs(phi_t),'linewidth',2);
xlabel('\tau');ylabel('|\phi(\tau)|');xlim([ta,tmx]);ylim([0,1.2]);grid on;
subplot(222);plot(t,Ek_t,'r',t,abs(phi_t.^2),'g',t,Ek_t+abs(phi_t.^2),'b','linewidth',2);
xlabel('\tau');ylabel('ENERGY');xlim([ta,tmx]);
legend('E_k','|\phi|^2','E_{total}',2);legend('boxoff');
subplot(223);plot(t,imag(Omega_t),'linewidth',2);hold on; 
plot([ta,tmx],[0,0],'r--');plot([0,tmx],[sqrt(3)/2,sqrt(3)/2],'r--');
xlabel('\tau');ylabel('\Omega_{IMAG}');xlim([ta,tmx]);
subplot(224);plot(t,real(Omega_t),'linewidth',2);hold on;
plot([0,tmx],[-0.5,-0.5],'r--');
xlabel('\tau');ylabel('\Omega_{REAL}');xlim([ta,tmx]);
print(gcf,'-dpng','nonlinearbeam-t_2');
