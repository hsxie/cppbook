% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2012-11-20 21:33
% half spectral method for mhd waves
% B0=1, rho0=1
% Alfven & fast & slow mode perfect match with theory
function mhd_waves
    close all; clear; clc;
    global k theta B0 rho0 va vs va2 vs2 sinth costh
    B0=1.0; rho0=1.0; va=2; vs=1; k=1.0; theta=pi/4;
    va2=va^2; vs2=vs^2;
    sinth=sin(theta); costh=cos(theta);
    y0=[0.1,-0.1,0.1,0.2,0.1,0.1,-0.1i];
    
    % three theory solutions, slow magnetosonic, alfven, fast magnetosonic
    w1=k*sqrt((va2+vs2)/2-sqrt((va2-vs2)^2+4*va2*vs2*sinth^2)/2);
    w2=k*sqrt(va2*costh^2);
    w3=k*sqrt((va2+vs2)/2+sqrt((va2-vs2)^2+4*va2*vs2*sinth^2)/2);
    
    [T,Y]=ode45(@push,0:0.02:1e2,y0);
    tt=T;
    drho=Y(:,1); dux=Y(:,2); duy=Y(:,3); duz=Y(:,4);
    dBx=Y(:,5); dBy=Y(:,6); dBz=Y(:,7);
    
    h=figure('unit','normalized','Position',[0.01 0.34 0.7 0.57]);
    set(gcf,'DefaultAxesFontSize',15);
    
    subplot(3,3,1:2); plot(tt,real(drho),tt,imag(drho),'LineWidth',2);
    xlabel('t'); ylabel('\delta\rho'); axis tight; grid on;
    title(['k=',num2str(k),', V_A=',num2str(va),', V_s=',num2str(vs),...
        ', \theta=',num2str(theta*180/pi),'^\circ']);
    subplot(334); plot(tt,real(dux),tt,imag(dux),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}ux'); axis tight; grid on;
    subplot(335); plot(tt,real(duy),tt,imag(duy),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}uy'); axis tight; grid on;
    subplot(336); plot(tt,real(duz),tt,imag(duz),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}uz'); axis tight; grid on;
    subplot(337); plot(tt,real(dBx),tt,imag(dBx),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}Bx'); axis tight; grid on;
    subplot(338); plot(tt,real(dBy),tt,imag(dBy),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}By'); axis tight; grid on;
    subplot(339); plot(tt,real(dBz),tt,imag(dBz),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}Bz'); axis tight; grid on;

    Lt=length(tt); % number of sampling
    dfs=2*pi/(tt(end)-tt(1));
    fs=0:dfs:dfs*(Lt-1);
    drho_ft=fft(real(drho))/Lt*2; % *2 ?? need check
    duy_ft=fft(real(duy))/Lt*2; % *2 ?? need check
    ifs=round(tt(end)*k*sqrt(vs2+va2)/4);
    subplot(333);
    plot(fs(1:ifs),abs(drho_ft(1:ifs)+abs(duy_ft(1:ifs))),'LineWidth',2);
%     title('simulation frequency');
    title(['\omega_{theory}=',num2str(w1),', ',...
        num2str(w2),', ',num2str(w3)]);
    ylabel('Amp');xlabel('\omega');
    xlim([0, fs(ifs)]); grid on;
    Amax=1.5*max(abs(drho_ft(1:ifs))); ylim([0,Amax]);
    hold on; plot([w1,w1],[0,Amax],'r--',[w2,w2],[0,Amax],'r--',...
        [w3,w3],[0,Amax],'r--','LineWidth',2);

end

function dy=push(t,y)
    global k B0 rho0 va2 vs2 sinth costh
    % y --> drho, dux, duy, duz, dBx, dBy, dBz
    drho=y(1);dux=y(2);duy=y(3);duz=y(4);dBx=y(5);dBy=y(6);dBz=y(7);
    dy=zeros(7,1);
    dy(1)=-1i*k*rho0*(dux*sinth+duz*costh);
    dy(2)=1i*k*va2*(costh*dBx/B0-sinth*dBz/B0)-1i*k*vs2*sinth*drho/rho0;
    dy(3)=1i*k*va2*costh*dBy/B0;
    dy(4)=-1i*k*vs2*costh*drho/rho0;
    dy(5)=1i*k*B0*dux*costh;
    dy(6)=1i*k*B0*duy*costh;
    dy(7)=-1i*k*B0*dux*sinth;
end
