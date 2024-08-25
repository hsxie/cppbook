% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2012-11-19 10:11
% half spectral method for em cold waves
% wce=1 unit
% match theory solution perfectly
function em_cold_waves
    close all; clear; clc;
    global k theta B0 qi qe mi me ni0 ne0 epsilon0 mu0 sinth costh ...
        wci wce wpi wpe c c2
    k=0.5; theta=pi/4; B0=1.0; qi=1; qe=-1; mi=4; me=1; ni0=1; ne0=1;
    epsilon0=4; mu0=1;
    c2=1/(epsilon0*mu0); c=sqrt(c2);
    sinth=sin(theta); costh=cos(theta);
    wci=qi*B0/mi; wce=qe*B0/me;
    wpi=sqrt(ni0*qi^2/epsilon0/mi); wpe=sqrt(ne0*qe^2/epsilon0/me);
    
    y0=[0.1,-0.1,0.2,0.1,0,0.1,0,-0.1,0.01,0.1,0,0];
    
    [T,Y]=ode45(@push,0:0.02:5e2,y0);
    tt=T;
    dEx=Y(:,1); dEy=Y(:,2); dEz=Y(:,3);
    dBx=Y(:,4); dBy=Y(:,5); dBz=Y(:,6);
    dvex=Y(:,7); dvey=Y(:,8); dvez=Y(:,9);
    dvix=Y(:,10); dviy=Y(:,11); dviz=Y(:,12);
    
    h=figure('unit','normalized','Position',[0.01 0.32 0.87 0.59]);
    set(gcf,'DefaultAxesFontSize',15);
    
    subplot(3,3,1:2); plot(tt,real(dEx),tt,imag(dEx),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}Ex'); axis tight; grid on;
    title(['k=',num2str(k),', \omega_{ce}=',num2str(wce),...
        ', \omega_{pe}=',num2str(wpe),', m_i/m_e=',num2str(mi/me),...
        ', \theta=',num2str(theta*180/pi),'^\circ']);
    subplot(334); plot(tt,real(dvix),tt,imag(dvix),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}u{ix}'); axis tight; grid on;
    subplot(335); plot(tt,real(dviy),tt,imag(dviy),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}u{iy}'); axis tight; grid on;
    subplot(336); plot(tt,real(dviz),tt,imag(dviz),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}u{iz}'); axis tight; grid on;
    subplot(337); plot(tt,real(dBx),tt,imag(dBx),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}Bx'); axis tight; grid on;
    subplot(338); plot(tt,real(dBy),tt,imag(dBy),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}By'); axis tight; grid on;
    subplot(339); plot(tt,real(dBz),tt,imag(dBz),'LineWidth',2);
    xlabel('t'); ylabel('\delta{}Bz'); axis tight; grid on;

    %% theory solutions
    
    wp=sqrt(wpe*wpe+wpi*wpi);

    kc=k*c;
    thetap=theta;

    % To speed up the program somewhat introduce these
    kc2=kc*kc;
    kc4=kc2*kc2;
    wci2=wci*wci;
    wp2=wp*wp;
    wextra=wp2+wci; % wce=1.0
    wextra2=wextra*wextra;
    cos2theta=cos(thetap)*cos(thetap);
    sin2theta=sin(thetap)*sin(thetap);

    % polynomial coefficients for EM model
    emc8=-(2*kc2+1+wci2+3*wp2);
    emc6=(kc4+(2*kc2+wp2)*(1+wci2+2*wp2)+wextra2);
    emc4=-(kc4*(1+wci2+wp2)+2*kc2*wextra2...
        +kc2*wp2*(1+wci2-wci)*(1+cos2theta)+wp2*wextra2);
    emc2=(kc4*(wp2*(1+wci2-wci)*cos2theta+wci*wextra)+kc2*wp2* ...
         wci*wextra.*(1+cos2theta));
    emc0=-kc4*wci2*wp2*cos2theta;

    empoly=[1, 0, emc8, 0, emc6, 0, emc4, 0, emc2, 0, emc0];
    wtmpem=roots(empoly); wem=sort(wtmpem);
    w1=wem(6);w2=wem(7);w3=wem(8);w4=wem(9);w5=wem(10);
  
    Lt=length(tt); % number of sampling
    dfs=2*pi/(tt(end)-tt(1));
    fs=0:dfs:dfs*(Lt-1);
    dEx_ft=fft(real(dEx))/Lt*2; % *2 ?? need check
    dEz_ft=fft(real(dEz))/Lt*2; % *2 ?? need check
    ifs=round(tt(end)*max(wem)/5);
    subplot(333);
    plot(fs(1:ifs),abs(dEx_ft(1:ifs))+abs(dEz_ft(1:ifs)),'LineWidth',2);
%     title('simulation frequency');
%     title(['\omega_{theory}=',num2str(w1),', ',num2str(w2),', ',num2str(w3)]);
    ylabel('Amp');xlabel('\omega');
    xlim([0, fs(ifs)]); grid on;    
    
    Amax=1.5*max(abs(dEx_ft(1:ifs))+abs(dEz_ft(1:ifs))); 
    subplot(333); ylim([0,Amax]);
    hold on; plot([w1,w1],[0,Amax],'r--',[w2,w2],[0,Amax],'r--',...
        [w3,w3],[0,Amax],'r--',[w4,w4],[0,Amax],'r--',...
        [w5,w5],[0,Amax],'r--','LineWidth',2);
    title(['\omega_{theory}=',num2str(w1),', ',num2str(w2),...
      ', ',num2str(w3),', ',num2str(w4),', ',num2str(w5)]);  
end

function dy=push(t,y)
    global k theta B0 qi qe mi me ni0 ne0 epsilon0 mu0 sinth costh ...
        wci wce wpi wpe c c2
    % y -> dEx, dEy, dEz, dBx, dBy, dBz, dvex, dvey, dvez, dvix, dviy, dviz
    dEx=y(1);dEy=y(2);dEz=y(3);dBx=y(4);dBy=y(5);dBz=y(6);
    dvex=y(7);dvey=y(8);dvez=y(9);dvix=y(10);dviy=y(11);dviz=y(12);
    dy=zeros(12,1);
    dJx=ni0*qi*dvix+ne0*qe*dvex;
    dJy=ni0*qi*dviy+ne0*qe*dvey;
    dJz=ni0*qi*dviz+ne0*qe*dvez;
    dy(1)=-1i*k*c2*dBy*costh-dJx/epsilon0;
    dy(2)=1i*k*c2*(dBx*costh-dBz*sinth)-dJy/epsilon0;
    dy(3)=1i*k*c2*dBy*sinth-dJz/epsilon0;
    dy(4)=1i*k*dEy*costh;
    dy(5)=-1i*k*(dEx*costh-dEz*sinth);
    dy(6)=-1i*k*dEy*sinth;
    dy(7)=(dEx+dvey*B0)*qe/me;
    dy(8)=(dEy-dvex*B0)*qe/me;
    dy(9)=dEz*qe/me;
    dy(10)=(dEx+dviy*B0)*qi/mi;
    dy(11)=(dEy-dvix*B0)*qi/mi;
    dy(12)=dEz*qi/mi;
end
