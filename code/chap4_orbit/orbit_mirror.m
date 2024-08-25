% Hua-sheng XIE, huashengxie@gmail.com, 2015-06-13 21:20
% Single particle orbit in dipole field, rewrite Orbit.m

close all; clear; clc;

global c q m Bfield rel Efield;

c0=2e-7; % mu0/2pi=2e-7
I=1e6; % current, unit: A
a=0.5; % radius of current coil, unit: m
L=4.0; % unit: m
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

figure('units','normalized','position',[0.01,0.07,0.55,0.75],...
    'DefaultAxesFontSize',15);

rel=0; % rel=1, relativistic else non-relativistic

% Parameters for dipole (Earth) and double dipole (Earth) B field
e = 1.602176565e-19; % Elementary charge (Coulomb)
m_pr = 1.672621777e-27; % Proton mass (kg)
m_el = 9.10938291e-31; % Electron mass (kg)
c = 299792458; % speed of light (m/s)

% Defintions for different E and B fields
Ezero=@(x,y,z,t)[0,0,0];
Efield=@(x,y,z,t)Ezero(x,y,z,t);

Bmirror=@(x,y,z,t)[fBx(x,y,z-L/2)+fBx(x,y,z+L/2),...
    fBy(x,y,z-L/2)+fBy(x,y,z+L/2),...
    fBz(x,y,z-L/2)+fBz(x,y,z+L/2)];

q=e; m=m_pr;
B0 = 3.07e-5; % Tesla
Bfield=@(x,y,z,t)Bmirror(x,y,z,t);

for ipl=1:2

if(ipl==1)
    axes('position',[0.09, 0.07, 0.4 ,0.8]); q=e; m=m_pr;
else
    axes('position',[0.59, 0.07, 0.4 ,0.8]); %q=-e; m=m_el;
    q=e; m=m_pr;
end

orbit=ipl;
if(orbit==1) % Trajectory of a proton with 10 MeV kinetic energy        
    K=1e2; % kinetic energy in eV
    x0=0; y0=0.8; z0=0;
    pitch_angle=25.0; % initial angle between velocity and mag.field (degrees)
    nT=100;
else
    K=1e3; % kinetic energy in eV
    x0=0; y0=2.8; z0=0;
    pitch_angle=5.0;
    nT=10;
end

K=K*e  ; % convert to Joule
% v=c/sqrt(1+(m_pr*c^2)/K); % replace m_pr with m_el for electron
v=c/sqrt(1+(m*c^2)/K); % replace m_pr with m_el for electron

if pitch_angle==90
    vz0=0 ; % avoids some numerical instabilities
else
    vz0=v*cos(pitch_angle*pi/180);
end
vy0=v*sin(pitch_angle*pi/180); vx0=0;

gamma=1.0/sqrt(1-(vx0^2+vy0^2+vz0^2)/c^2);
T=2*pi*gamma*m/(abs(q)*sqrt(sum(Bfield(x0,y0,z0,0).^2)));
dt=T/50;
tend=nT*T; % end time

yy0=[x0,y0,z0,vx0,vy0,vz0];
[t,y]=ode45('SolveNewtonLorenz',0:dt:tend,yy0);

%% Plotting
cla;
plot3(x0,y0,z0,'ro','Linewidth',2); hold on; 
% [X,Y,Z]=sphere(20); surf(X,Y,Z); hold on; grid on;
plot3(y(:,1),y(:,2),y(:,3),'Linewidth',2);
xlabel('x'); ylabel('y'); zlabel('z');
amax=2*sqrt(x0^2+y0^2+z0^2); amin=-amax;

hold on; ang=0:pi/50:2*pi;
plot3(a.*cos(ang),a.*sin(ang),0.*cos(ang)+L/2,'r','LineWidth',5);
hold on;
plot3(a.*cos(ang),a.*sin(ang),0.*cos(ang)-L/2,'r','LineWidth',5);
axis equal; box on;
% axis off;

title([num2str(K/e,'%2.1e'),'eV, m=',num2str(m/m_el,3),'m_e',10,'pitch angle ',...
    num2str(pitch_angle),'^o']);

end

set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng','orbit_dipole.png');
% print(gcf,'-dpdf','-painters','orbit_dipole.pdf');
% saveas(gcf,'orbit_dipole.fig','fig');
