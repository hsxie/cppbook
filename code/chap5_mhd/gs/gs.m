% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-03-27 15:10
% Solve Grad-Shafranov equation for Tokamak equilibrium, with fixed
% boundary
% The treatment of the boundary grid is still too mess. Any better one?
close all; clear; clc;
runtime=cputime;

% keep ap>=1 and aI>=0.5
% R0=1.0; a=0.3; k=1.5; d=0.5; ap=2.0; aI=2.5;
% R0=3.0; a=1.0; k=1.2; d=0.5; ap=2.0; aI=3.0;
% R0=3.0; a=1.0; k=1.2; d=0.5; ap=1.0; aI=1.0;
R0=3.0; a=1.0; k=1.0; d=0.0; ap=2.0; aI=2.0; % circle
% R0=3.0; a=1.0; k=1.2; d=0.5; ap=2.0; aI=0.6;

np=2; nI=2;

mu0=4.0e-8*pi; p0=0.01e6; I0=1e6;

nR=51; nZ=51; Rlow=R0-a; Rup=R0+a; Zlow=0; Zup=k*a;

rr=linspace(Rlow,Rup,nR); zz=linspace(Zlow,Zup,nZ);
dR=rr(2)-rr(1); dZ=zz(2)-zz(1);

[R,Z]=meshgrid(rr,zz); dR2=dR*dR; dZ2=dZ*dZ;

%%
psi_n=0:0.02:1; p=p0.*(1-psi_n.^np).^ap; I=I0.*(1-psi_n.^nI).^aI;

hf=figure('unit','normalized','Position',[0.01 0.1 0.65 0.65],...
            'DefaultAxesFontSize',15);
subplot(231); plot(psi_n,p,'LineWidth',2); xlabel('\psi_n'); ylabel('p'); 
title(['p=p_0(1-\psi_n^{',num2str(np),'})^{',num2str(ap),'}']);
subplot(232); plot(psi_n,I,'LineWidth',2);  xlabel('\psi_n'); ylabel('I');
title(['I=I_0(1-\psi_n^{',num2str(nI),'})^{',num2str(aI),'}']);

% plot grid
subplot(233);
plot(repmat(rr,length(zz),1),repmat(zz.',1,length(rr)),'b',...
    repmat(rr.',1,length(zz)),repmat(zz,length(rr),1),'b');
 
% cal jend
the=0:pi/200:pi;
Rb=R0+a.*cos(the+d.*sin(the)); % boundary
Zb=k*a.*sin(the);
Rbg=rr; % R boundary on grid
Zbg=interp1(Rb,Zb,Rbg); % interp
jend=floor((Zbg-Zlow)./dZ)+1;
Zbg=Zlow+(jend-1).*dZ; % get the Z boundary on grid

hold on; plot(Rb,Zb,Rbg,Zbg,'r+','LineWidth',2);

% cal iend
the=0:pi/100:pi/2;
Rb=R0+a.*cos(the+d.*sin(the)); % boundary
Zb=k*a.*sin(the);
Zbg=zz; % Z boundary on grid
Rbg=interp1(Zb,Rb,Zbg); % interp
iend=floor((Rbg-Rlow)./dR)+1;
Rbg=Rlow+(iend-1).*dR; % get the Z boundary on grid
hold on;plot(Rbg,Zbg,'gx','LineWidth',2);

% cal istart
the=(pi/2-eps):pi/100:(pi+eps);
Rb=R0+a.*cos(the+d.*sin(the)); % boundary
Zb=k*a.*sin(the);
Zbg=zz; % Z boundary on grid
Rbg=interp1(Zb,Rb,Zbg,'cubic'); % interp
istart=ceil((Rbg-Rlow)./dR)+1;
istart(find(istart<=1))=2;
Rbg=Rlow+(istart-1).*dR; % get the Z boundary on grid
hold on; plot(Rbg,Zbg,'gx','LineWidth',2);
xlabel('R'); ylabel('Z'); axis equal; axis tight; 
title('grid');

%% main
psib=0.0; % set psib=0
psi0=-10; % initial guess
psi0n=psi0+1;

ind=find(((R-R0).^2+(Z./k).^2)<a^2/4); % initial
psi=0.*R+psib;
psi(ind)=psi0.*exp(-(((R(ind)-R0)./a).^2+(Z(ind)./(k*a)).^2)*5);
% S=0.*psi;
itr=0;
while(itr<=40 && abs(psi0n-psi0)>1e-3) % iteration, for psi_n 
    psi0n=psi0;
    nt=1; % iterative times
    psim1=max(max(abs(psi)));
    psim2=psim1+1;
    while(nt<=1000 && abs(psim1-psim2)/psim1>1e-3) % iteration, G-S eq
%     while(nt<=400)
        psim1=psim2;
        psi_bar=(psi-psi0)./(psib-psi0);
        S=(mu0.*R.^2.*p0.*ap.*np.*psi_bar.^(np-1).*(1-...
            psi_bar.^np).^(ap-1)+(mu0/2/pi*I0)^2.*aI.*nI.*psi_bar.^(nI-...
            1)*(1-psi_bar.^nI).^(2*aI-1))/(psib-psi0);

        for i=2:nR-1
            for j=2:jend(i)-1
                if((i>=istart(j))&&(i<=iend(j)))
                    psi(j,i)=((psi(j,i-1)+psi(j,i+1))/dR2+...
                        (psi(j-1,i)+psi(j+1,i))/dZ2...
                        +(psi(j,i-1)-psi(j,i+1))/(2*dR)/R(j,i)-...
                        S(j,i))/(2.0/dR2+2.0/dZ2);
                end
            end
            psi(1,i)=psi(2,i);
        end
        nt=nt+1;
        psim2=max(max(abs(psi)));
    end
%     psi0=max(max((psi)));
    psi0=min(min((psi)));
    itr=itr+1;
end
psi=psi-psi0; % reset psi to keep the on axis psi be zero
psib=psib-psi0;
indm=find(psi==0);
Rm=min(R(indm)); % find axis
Zm=min(Z(indm));

%% plot
v=psi0:(psib-psi0)/10:psib;
subplot(234); contour(R,Z,psi,10,'LineWidth',2); 
hold on; plot(Rm,Zm,'rx','LineWidth',2,'MarkerSize',5);
xlabel('R'); ylabel('Z'); axis equal; axis tight;
title(['R0=',num2str(R0),', a=',num2str(a),', k=',...
    num2str(k),', d=',num2str(d)]);

subplot(235); surf(R,Z,psi); 
xlabel('R'); ylabel('Z'); zlabel('\psi');
subplot(236); pcolor(R,Z,psi); shading interp;
xlabel('R'); ylabel('Z'); axis equal; axis tight;

runtime=cputime-runtime; title(['Run time: ',num2str(runtime),'s']);

print('-dpng',['gs_plot_R0=',num2str(R0),',a=',num2str(a),...
    ',k=',num2str(k),',d=',num2str(d),',np=',num2str(np),...
    ',nI=',num2str(nI),',ap=',num2str(ap),',aI=',num2str(aI),'.png']);

