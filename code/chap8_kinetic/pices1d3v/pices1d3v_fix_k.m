close all;clear all;clc;

kwc=2.5; % k*rho_c

wp=1; qm=-1; wc=1/2.5; theta = pi/2; B0=wc/qm; Bx0=B0*cos(theta); 
By0=B0*sin(theta);

k=kwc*wc; % k*lambda_D

% parameters
L=2*pi/k; dt=.05; nt=4000; ntout=nt/2; ng=32; np=10000;
vb=1.0; xp1=5.0e-1; vp1=0.0;
vt=1; % note: the normalization sqrt(2) will be found in randn()
q=wp^2/(qm*np/L); rho_back=-q*np/L; dx=L/ng;

% initial loading for the 2 Stream instability
xp=linspace(0,L,np)';
% vp=vt*randn(np,1)+(1-2*mod([1:np]',2)).*vb;
vpx=vt*randn(np,1); % randn is {exp[-(x-mu)^2/(2*sgm^2)]}/[sgm*sqrt(2*pi)]
vpy=vt*randn(np,1);
vpz=vt*randn(np,1);

% Perturbation
% vpx=vpx+vp1*cos(k*xp);
xp=xp+xp1*cos(k*xp);
p=1:np;p=[p p];

% Main computational cycle
h = figure('Unit','Normalized','position',...
    [0.02 0.1 0.6 0.7],'DefaultAxesFontSize',15);
for it=1:nt
    % apply periodic bc on the particle positions
    xp=xp./L+10.0; xp=L.*(xp-floor(xp));
    
    % diagnosing
    if(mod(it,ntout)==0)
        subplot(221);plot(xp,vpx,'b.','Markersize',2);
        axis([0,L,-3*(abs(vt)+abs(vb)),3*(abs(vt)+abs(vb))]);
        title(['(a) phase space, t=',num2str(it*dt),', np=',num2str(np)]);
        xlabel('xp');ylabel('vpx'); pause(0.2);
%         print(gcf, '-dpng', ['vpx-x,t=',num2str(it*dt),'.png']);
    end
    
    % update xp
    xp=xp+vpx*dt;
    
    % projection p->g
    g1=floor(xp/dx-.5)+1;g=[g1;g1+1];
    fraz1=1-abs(xp/dx-g1+.5);
    fraz=[fraz1;1-fraz1];
    
    % apply bc on the projection
    out=(g<1);g(out)=g(out)+ng;
    out=(g>ng);g(out)=g(out)-ng;
    mat=sparse(p,g,fraz,np,ng);
    rho=full((q/dx)*sum(mat))'+rho_back;
    
    % computing fields, dE/dx
    Eg=zeros(ng,1);
    for j=1:ng-1
        Eg(j+1)=Eg(j)+(rho(j)+rho(j+1))*dx/2;
    end
    Eg(1)=Eg(ng)+rho(ng)*dx;
    Eg=Eg-mean(Eg);
    
    % projection q->p and update of vp
    vpx=vpx+qm*(mat*Eg+vpy*0-vpz*By0)*dt;
    vpy=vpy+qm*(vpz*Bx0-vpx*0)*dt;
    vpz=vpz+qm*(vpx*By0-vpy*Bx0)*dt;
    
    EEk(it)=0.5*abs(q)*sum(vpx.^2+vpy.^2+vpz.^2); % kinetic energy
    EEf(it)=0.5*sum(Eg.^2)*dx; % potential energy
    Et(it)=Eg(floor(ng/2));
    t(it)=it*dt;
end
%%
subplot(222);plot(t,EEk,t,EEf,t,EEk+EEf,'r:','LineWidth',2);
title(['(b) \omega_c/\omega_p=',num2str(wc),', nt=',num2str(nt)]);
xlabel('t'); ylabel('Energy');legend('E_k','E_e','E_{tot}',4);
legend('boxoff');xlim([0,max(t)]);

% subplot(122);

% % Find the corresponding indexes of the extreme max values
% lndE=log(sqrt(real((EEf(1:nt))))); % EEf=E^2=[exp(gam*t)]^2=exp(2*gam*t)
% it0=floor(nt*1/20); it1=floor(nt*17/20);
% yy=lndE(it0:it1);
% extrMaxIndex = find(diff(sign(diff(yy)))==-2)+1;
% t1=t(it0+extrMaxIndex(1));t2=t(it0+extrMaxIndex(end));
% y1=yy(extrMaxIndex(1));y2=yy(extrMaxIndex(end));
% plot(t,lndE,[t1,t2],[y1,y2],'r*--','LineWidth',2);
% omega=pi/((t2-t1)/(length(extrMaxIndex)-1));
% gammas=(real(y2)-real(y1))/(t2-t1);
% title(['(b) \omega^S=',num2str(omega),', \gamma^S=',num2str(gammas)]);
% axis tight;

% subplot(222);plot(t,EEf,'LineWidth',2);

nt=length(t); y=Et; yf=fft(y); tf=1/dt.*linspace(0,1,nt); tf=tf/wc*2*pi; %
subplot(223);
plot(t,y,'LineWidth',2); title(['(c) Et(ng/2)',', dt=',num2str(dt)]);
xlabel('t\omega_p^{-1}');xlim([0,max(t)]);
subplot(224); plot(tf,abs(yf),'LineWidth',2);
for jp=1:10
    hold on; plot([jp,jp],[0,max(abs(yf))],'g--','LineWidth',2);
end
xlabel('\omega/\omega_c'); 
title(['(d) spectral, k\rho_c=',num2str(kwc),...
    ', \Delta\omega=',num2str(tf(2)-tf(1))]);
% xlim([0,tf(floor(nt/5))]);
xlim([0,10]); ylim([0,max(abs(yf))]);

print(gcf,'-dpng',['pices1d3v_np',num2str(np),'_wc',num2str(wc),...
    '_k',num2str(kwc),'_dt',num2str(dt),'_nt',num2str(nt),'.png']);


