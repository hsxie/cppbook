% Hua-sheng XIE, huashengxie@gmail.com, 2015-06-14 13:33
% Plot dipole field from single current coil

close all; clear; clc;

c0=2e-7; % mu0/2pi=2e-7
I=1e6; % current, unit: A
a=0.5; % radius of current coil, unit: m
eps=1e-4;

% Set field for single coil
rho=@(x,y,z)sqrt(x.^2+y.^2);
theta=@(x,y,z)atan2(y,x);
k2_0=@(x,y,z)4*a.*rho(x,y,z)./((a+rho(x,y,z)).^2+z.^2);
k2=@(x,y,z)k2_0(x,y,z)-eps*(k2_0(x,y,z)==1);
% K=@(x,y,z)KK(k2(x,y,z)); % [K,E]=ellipke(x);
% E=@(x,y,z)EE(k2(x,y,z));
K=@(x,y,z)pi/2.*(1+(1/2)^2*k2(x,y,z)+(1/2*3/4)^2*k2(x,y,z).^2+(1/2*3/4*5/6)^2*k2(x,y,z).^3);
E=@(x,y,z)pi/2.*(1-(1/2)^2*k2(x,y,z)+(1/2*3/4)^2*k2(x,y,z).^2/3-(1/2*3/4*5/6)^2*k2(x,y,z).^3/5);
a1=@(x,y,z)sqrt((a+rho(x,y,z)).^2+z.^2);
b1_0=@(x,y,z)(a-rho(x,y,z)).^2+z.^2;
b1=@(x,y,z)b1_0(x,y,z)+eps*(b1_0(x,y,z)==0);
b2=@(x,y,z)a^2+rho(x,y,z).^2+z.^2;
b3=@(x,y,z)a^2-rho(x,y,z).^2-z.^2;
Brho=@(x,y,z)c0*I.*z./((rho(x,y,z)+eps*(rho(x,y,z)==0))...
    .*a1(x,y,z)).*(b2(x,y,z)./b1(x,y,z).*E(x,y,z)-K(x,y,z));
Bscz=@(x,y,z)c0*I./a1(x,y,z).*(b3(x,y,z)./b1(x,y,z).*E(x,y,z)+K(x,y,z));
Btheta=@(x,y,z)0;

Bscx=@(x,y,z)Brho(x,y,z).*cos(theta(x,y,z))-Btheta(x,y,z).*sin(theta(x,y,z));
Bscy=@(x,y,z)Brho(x,y,z).*sin(theta(x,y,z))+Btheta(x,y,z).*cos(theta(x,y,z));

fBx=@(x,y,z)Bscx(x,y,z);
fBy=@(x,y,z)Bscy(x,y,z);
fBz=@(x,y,z)Bscz(x,y,z);
fB=@(x,y,z)sqrt(fBx(x,y,z).^2+fBy(x,y,z).^2+fBz(x,y,z).^2);

figure('units','normalized','position',[0.02,0.1,0.6,0.5],...
    'DefaultAxesFontSize',15);

ri=[1.1, 0.0, 0.0, 32;
    1.2, 0.0, 0.0, 63;
    1.5, 0.0, 0.0, 150;
    2.1, 0.0, 0.0, 300;
    3.0, 0.0, 0.0, 450;
    0.1, 0.0, 0.0, 200;
    0.2, 0.0, 0.0, 250;
    ];
ri=[ri;-ri];

xi=a*ri(:,1); yi=a*ri(:,2); zi=a*ri(:,3); nti=abs(ri(:,4)); 

npl=length(xi);
for theta=0:30:150;
    
    for ipl=1:npl

        x0=xi(ipl)*cos(theta*pi/180); y0=xi(ipl)*sin(theta*pi/180); z0=zi(ipl);
        dl=0.01/2; x=[]; y=[]; z=[]; nt=nti(ipl);
        x(1)=x0; y(1)=y0; z(1)=z0;
        for pm=[-1,1]
            for it=1:nt % need update to Boris/RK-4 ...
                x(it+1)=x(it)+dl*fBx(x(it),y(it),z(it))/fB(x(it),y(it),z(it))*pm;
                y(it+1)=y(it)+dl*fBy(x(it),y(it),z(it))/fB(x(it),y(it),z(it))*pm;
                z(it+1)=z(it)+dl*fBz(x(it),y(it),z(it))/fB(x(it),y(it),z(it))*pm;
            end
            subplot(121); plot3(x,y,z,'Linewidth',2); hold on; box on;
            if(theta==0)
                subplot(122); plot(x,z,'Linewidth',2); hold on; box on;
            end
        end
    end
end
%%
subplot(121); xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
hold on; ang=0:pi/50:2*pi;
plot3(a.*cos(ang),a.*sin(ang),0.*cos(ang),'ro','LineWidth',5);

subplot(122); xlabel('x'); ylabel('z'); axis equal; axis tight;
plot([-a,a],[0,0],'ro'); hold on; plot([-a,a],[0,0],'rx'); hold on;

set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng','B_coilloop.png');
print(gcf,'-dpdf','-painters','B_coilloop.pdf');
saveas(gcf,'B_coilloop.fig','fig');
