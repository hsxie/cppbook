function siglecoil_field_rhothetaz_fun2

close all;clear;clc;

c0=2e-7; % mu0/2pi=2e-7
I=1e6; % current, unit: A
a=0.2; L=2.0;
eps=1e-4;

% Set field for single coil
rho=@(x,y,z)sqrt(x.^2+y.^2);
theta=@(x,y,z)atan2(y,x);
k2_0=@(x,y,z)4*a.*rho(x,y,z)./((a+rho(x,y,z)).^2+z.^2);
k2=@(x,y,z)k2_0(x,y,z)-eps*(k2_0(x,y,z)==1);
K=@(x,y,z)KK(k2(x,y,z));
E=@(x,y,z)EE(k2(x,y,z));
a1=@(x,y,z)sqrt((a+rho(x,y,z)).^2+z.^2);
b1_0=@(x,y,z)(a-rho(x,y,z)).^2+z.^2;
b1=@(x,y,z)b1_0(x,y,z)+eps*(b1_0(x,y,z)==0);
b2=@(x,y,z)a^2+rho(x,y,z).^2+z.^2;
b3=@(x,y,z)a^2-rho(x,y,z).^2-z.^2;
Brho=@(x,y,z)c0*I.*z./((rho(x,y,z)+eps*(rho(x,y,z)==0)).*a1(x,y,z)).*(b2(x,y,z)./b1(x,y,z).*E(x,y,z)-K(x,y,z));
Bscz=@(x,y,z)c0*I.*z./a1(x,y,z).*(b3(x,y,z)./b1(x,y,z).*E(x,y,z)+K(x,y,z));
Btheta=@(x,y,z)0;

Bscx=@(x,y,z)Brho(x,y,z).*cos(theta(x,y,z))-Btheta(x,y,z).*sin(theta(x,y,z));
Bscy=@(x,y,z)Brho(x,y,z).*sin(theta(x,y,z))+Btheta(x,y,z).*cos(theta(x,y,z));

Bsc=@(x,y,z)[Bscx(x,y,z),Bscy(x,y,z),Bscz(x,y,z)];

% Plot
[X,Y,Z]=meshgrid(-2*a:a/20:2*a,-2*a:a/20:2*a,-1.5*L:L/20:1.5*L);
% [Bx,By,Bz]=Bsc(X,Y,Z);
% Bx=Bscx(X,Y,Z);
% By=Bscy(X,Y,Z);
% Bz=Bscz(X,Y,Z);
Bx=Bscx(X,Y,Z+L/2)+Bscx(X,Y,Z-L/2);
By=Bscy(X,Y,Z+L/2)+Bscy(X,Y,Z-L/2);
Bz=Bscz(X,Y,Z+L/2)-Bscz(X,Y,Z-L/2);


% [rho,theta,z]=meshgrid(0.01*a:0.25*a:0.9*a,0:pi/6:2*pi,0.01);
[rho,theta,z]=meshgrid(0.01*a:0.02*a:0.1*a,0:pi/6:2*pi,-0.45*L:0.1*L:0.1*L);
sx=rho.*cos(theta); sy=rho.*sin(theta); sz=z;

% subplot(211);streamline(X,Y,Z,Bx,By,Bz,sx,sy,sz,[0.1,5000]);
% axis equal;zlim([-L,L]);view(3);camroll(270);
% xlabel('x');ylabel('y');zlabel('z');
% hold on; ang=0:pi/50:2*pi;
% plot3(a.*cos(ang),a.*sin(ang),0.*cos(ang)+L/2,'r','LineWidth',5);
% hold on;
% plot3(a.*cos(ang),a.*sin(ang),0.*cos(ang)-L/2,'r','LineWidth',5);

% subplot(212);
streamline(X,Y,Z,Bx,By,Bz,sx,sy,sz,[0.1,5000]);
axis equal;zlim([-L,L]);view(90,0);camroll(270);
xlabel('x');ylabel('y');zlabel('z');
hold on; ang=0:pi/50:2*pi;
plot3(a.*cos(ang),a.*sin(ang),0.*cos(ang)+L/2,'r','LineWidth',5);
hold on;
plot3(a.*cos(ang),a.*sin(ang),0.*cos(ang)-L/2,'r','LineWidth',5);

end

%% Ellipse function
function K=KK(x)
    [K,E]=ellipke(x);
end
function E=EE(x)
    [K,E]=ellipke(x);
end