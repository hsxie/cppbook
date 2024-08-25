% Orbit.m, Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2012-01-04 21:22
% Single particle motion in magnetic fields and guiding-center trajectory
% E(x,y,z,t), B(x,y,z,t), with and without relativistic effect.
% Suggest Refs:
% [1] Ozturk, M. K. Trajectories of charged particles trapped in Earth's 
%   magnetic field arXiv, 2011, http://arxiv.org/abs/1112.3487.
% [2] http://ptsg.eecs.berkeley.edu/, S_PARMOS code.
% Note: The Guiding center (GC) version is not perfect yet, then not here.

% 2012-03-19 19:19 and 2012-12-28 10:10, Bugs in tokamak part are fixed.
%%
close all; clear; clc; 

global c emtype q m Bfield rel Efield;
global B0 Bn d R0 Re E0;

% select field: emtype=1, dipole (Earth); emtype=2, double dipole (Earth);
%          emtype=3, parabolic; emtype=4, Harris sheet; emtype=5, Tokamak;
%          emtype=6, E=[0,E0,0], B=[0,0,B0], for EXB drift.
%          emtype=7, E=[E0*cos(x),0,0], B=[0,0,B0], Non-uniform E.
%          else, uniform.
emtype=5;

rel=0; % rel=1, relativistic else non-relativistic

% Parameters for dipole (Earth) and double dipole (Earth) B field
e = 1.602176565e-19; % Elementary charge (Coulomb)
m_pr = 1.672621777e-27; % Proton mass (kg)
m_el = 9.10938291e-31; % Electron mass (kg)
c = 299792458; % speed of light (m/s)
Re = 6378137; % meter (Earth radius)

%% Defintions for different E and B fields
Ezero=@(x,y,z,t)[0,0,0];
Eonx=@(x,y,z,t)[E0,0,0];
Eony=@(x,y,z,t)[0,E0,0];
Eonz=@(x,y,z,t)[0,0,E0];
Ecosx=@(x,y,z,t)[E0*cos(x),0,0];

Efield=@(x,y,z,t)Ezero(x,y,z,t);

Buniform=@(x,y,z,t)[0,0,B0];

Bdipole=@(x,y,z,t)[3*x*z*(-B0*Re^3/power(x*x+y*y+z*z,5.0/2.0)),...
    3*y*z*(-B0*Re^3/power(x*x+y*y+z*z,5.0/2.0)),...
    (2*z*z-x*x-y*y)*(-B0*Re^3/power(x*x+y*y+z*z,5.0/2.0))];

Bdoubledipole=@(x,y,z,t)Bdipole(x,y,z,t)+Bdipole(x-20*Re,y,z);

Bparabolic=@(x,y,z,t)[(heaviside(1-abs(z))*z+heaviside(abs(z)-1)*sign(z))*B0/d,0,Bn];

Bharrissheet=@(x,y,z,t)[B0*tanh(z/d), 0, Bn];

rr=@(x,y,z,t)sqrt((sqrt(x^2+y^2)-R0)^2+z^2)+1e-10;
qq=@(x,y,z,t)1.0+rr(x,y,z,t)^2; % q profile
% Bp=@(x,y,z,t)B0/sqrt(1+(qq(x,y,z,t)*R0)^2/rr(x,y,z,t)^2); % not correct
% Bt=@(x,y,z,t)B0/sqrt(1+rr(x,y,z,t)^2/(qq(x,y,z,t)*R0)^2);
Bt=@(x,y,z,t)B0*R0/sqrt(x^2+y^2);
Bp=@(x,y,z,t)Bt(x,y,z,t)*rr(x,y,z,t)/(qq(x,y,z,t)*R0);
Btkmkx=@(x,y,z,t)(-Bt(x,y,z,t)*(y/sqrt(x^2+y^2))-Bp(x,y,z,t)*(z/rr(x,y,z,t))*(x/sqrt(x^2+y^2)));
Btkmky=@(x,y,z,t)(Bt(x,y,z,t)*(x/sqrt(x^2+y^2))-Bp(x,y,z,t)*(z/rr(x,y,z,t))*(y/sqrt(x^2+y^2)));
% Btkmkz=@(x,y,z,t)(Bp(x,y,z,t)*(sqrt(rr(x,y,z,t)^2-z^2)/rr(x,y,z,t)));% not correct, only +
Btkmkz=@(x,y,z,t)(Bp(x,y,z,t)*((sqrt(x^2+y^2)-R0)/rr(x,y,z,t)));
BBtkmk=@(x,y,z,t)sqrt(Btkmkx(x,y,z,t)^2+Btkmky(x,y,z,t)^2+Btkmkz(x,y,z,t)^2);
Btokamak=@(x,y,z,t)[Btkmkx(x,y,z,t),Btkmky(x,y,z,t),Btkmkz(x,y,z,t)];

if(emtype==1) % Dipole
    q=e;m=m_pr;
    B0 = 3.07e-5; % Tesla
    Bfield=@(x,y,z,t)Bdipole(x,y,z,t);
elseif(emtype==2) % Double dipole
    m=m_pr; q=e;
    B0 = 3.07e-5; % Tesla
    Bfield=@(x,y,z,t)Bdoubledipole(x,y,z,t);
elseif(emtype==3) % Parabolic
    B0=10.0; Bn=1.0; d=0.2; m=5; q=1;
    Bfield=@(x,y,z,t)Bparabolic(x,y,z,t);
elseif(emtype==4) % Harris Sheet
    B0=10.0; Bn=1.0; d=1; m=5; q=1;
    Bfield=@(x,y,z,t)Bharrissheet(x,y,z,t);
elseif(emtype==5) % Tokamak
    B0=1.0; R0=3.0; m=m_pr; q=e;
    Bfield=@(x,y,z,t)Btokamak(x,y,z,t);
elseif(emtype==6) % EXB
    B0=1.0; m=1; q=1; E0=1.0;
    Bfield=@(x,y,z,t)Buniform(x,y,z,t);
    Efield=@(x,y,z,t)Eony(x,y,z,t);
elseif(emtype==7) % non-uniform E
    B0=1.0; m=1; q=1; E0=1.0;
    Bfield=@(x,y,z,t)Buniform(x,y,z,t);
    Efield=@(x,y,z,t)Ecosx(x,y,z,t);
else % Uniform
    B0=1.0; m=1; q=1;
    Bfield=@(x,y,z,t)Buniform(x,y,z,t);
end

%% Initial and Calculating
if(emtype==1||emtype==2) % Dipole
    orbit=1;
    if(orbit==1) % Trajectory of a proton with 10 MeV kinetic energy        
        K=1e7; % kinetic energy in eV
        x0=4*Re; y0=0*Re; z0=0*Re;
        pitch_angle=30.0; % initial angle between velocity and mag.field (degrees)
    else
        K=2e5; % kinetic energy in eV
        x0=7*Re; y0=7*Re; z0=0*Re;
        pitch_angle=60.0;
    end
    
    K=K*e  ; % convert to Joule
    v=c/sqrt(1+(m_pr*c^2)/K); % replace m_pr with m_el for electron
    
    if pitch_angle==90
        vz0=0 ; % avoids some numerical instabilities
    else
        vz0=v*cos(pitch_angle*pi/180);
    end
    vy0=v*sin(pitch_angle*pi/180); vx0=0;
    
    gamma=1.0/sqrt(1-(vx0^2+vy0^2+vz0^2)/c^2);
    T=2*pi*gamma*m/(abs(q)*sqrt(sum(Bfield(x0,y0,z0,0).^2)));
    dt=T/10;
    tend=1000*T; % end time
elseif(emtype==3)
    orbit=1;
    if(orbit==1) % ring orbit
        pa = 30*pi/180;
        v = 2;
        x0=0.0; y0=0; z0=0.0;
        phi=-60*pi/180;
        vx0=v*sin(pa)*cos(phi); vy0 = v*sin(pa)*sin(phi); vz0 = v*cos(pa);
        tend=100;
    elseif(orbit==2) % speiser orbit
        x0=5.0; y0=-5; z0=1.0;
        vx0=-1; vy0=1; vz0=0;
        tend=30;
    elseif(orbit==3) % cucumber orbit
        x0=5.0; y0=0; z0=0.9;
        vx0=0; vy0=0.2; vz0=0.5;
        tend=300;
    else % test
        x0=1.0; y0=0; z0=0.1;
        vx0=0.1; vy0=-0.2; vz0=0.0;
        tend=300;
    end
    dt=0.01;
elseif(emtype==4)
    x0=5.0; y0=-5; z0=1.0;
    vx0=-1; vy0=1; vz0=0;
    tend=70;
    dt=0.01;
elseif(emtype==5)
    orbit=2;
    if(orbit==1) % Trajectory of a proton with 1 MeV kinetic energy        
        K=1e6; % kinetic energy in eV
        x0=(1.0+0.2)*R0; y0=0*R0; z0=0*R0;
        pitch_angle=30.0; % initial angle between velocity and mag.field (degrees)
    elseif(orbit==2)
        K=1e5; % kinetic energy in eV
        x0=(1.0+0.4)*R0; y0=0*R0; z0=0*R0;
        pitch_angle=0.0;
    else
        K=2e5; % kinetic energy in eV
        x0=1.2*R0; y0=0*R0; z0=0*R0;
        pitch_angle=30.0;
    end    
    K=K*e  ; % convert to Joule
    v=c/sqrt(1+(m_pr*c^2)/K); % replace m_pr with m_el for electron
    vpara0=v*cos(pitch_angle*pi/180);
    vperp0=v*sin(pitch_angle*pi/180);
    Btkmkx0=Btkmkx(x0,y0,z0,0);Btkmky0=Btkmky(x0,y0,z0,0);
    Btkmkz0=Btkmkz(x0,y0,z0,0);BBtkmk0=BBtkmk(x0,y0,z0,0);
    
    vx0=(vpara0*Btkmkx0+vperp0*Btkmkx0*Btkmkz0/sqrt(Btkmkx0^2+Btkmky0^2))/BBtkmk0;
    vy0=(vpara0*Btkmky0+vperp0*Btkmky0*Btkmkz0/sqrt(Btkmkx0^2+Btkmky0^2))/BBtkmk0;
    vz0=(vpara0*Btkmkz0-vperp0*sqrt(Btkmkx0^2+Btkmky0^2))/BBtkmk0;
    
    T=2*pi*m/(abs(q)*B0);
    dt=T/16;
    tend=400*T;
else
    x0=1.0; y0=0.0; z0=0.0;
    vx0=0; vy0=0.1; vz0=0.01;
    T=2*pi*m/(abs(q)*B0);
    dt=T/50;
    tend=6*T;
end

yy0=[x0,y0,z0,vx0,vy0,vz0];
options = odeset('RelTol',1e-4); % 17-04-07 13:38
[t,y]=ode45('SolveNewtonLorenz',0:dt:tend,yy0,options);

%% Plotting
plot3(x0,y0,z0,'o'); hold on; 
plot3(y(:,1),y(:,2),y(:,3),'r');
xlabel('x');ylabel('y');zlabel('z');
if(emtype==1)
    cla;
    [X,Y,Z]=sphere(20);surf(X,Y,Z);hold on; grid on;
    plot3(y(:,1)./Re,y(:,2)./Re,y(:,3)./Re); title('Dipole (Earth)');
    xlabel('x[R_e]');ylabel('y[R_e]');zlabel('z[R_e]');
    amax=2*sqrt(x0^2+y0^2+z0^2)/Re;amin=-amax;
    axis equal;
%     axis([amin,amax,amin,amax,amin,amax]);
elseif(emtype==2)
    cla;
    [X,Y,Z]=sphere(20);surf(X,Y,Z);hold on;
    plot3(y(:,1)./Re,y(:,2)./Re,y(:,3)./Re); title('Double dipole (Earth)');
    xlabel('x[R_e]');ylabel('y[R_e]');zlabel('z[R_e]');
    amax=2*sqrt(x0^2+y0^2+z0^2)/Re;amin=-amax;
    axis equal;
%     axis([amin,amax,amin,amax,amin,amax]);
elseif(emtype==3)
    title('Parabolic');
elseif(emtype==4)
    title('Harris Sheet');
elseif(emtype==5)
    hold on;
    v0=sqrt(vx0^2+vy0^2+vz0^2);
    quiver3(x0,y0,z0,vx0/v0,vy0/v0,vz0/v0,1.5);
    title(['Tokamak, R_0=',num2str(R0),'m, E=',num2str(K/e/1e3),'keV']);
    [u,v]=meshgrid(0:2*pi/50:2*pi,0:2*pi/50:1.5*pi);
    X=(R0+0.5.*R0.*cos(u)).*cos(v); Y=(R0+0.5.*R0.*cos(u)).*sin(v);
    Z=0.5.*R0.*sin(u); mesh(X,Y,Z);axis equal;hidden off;colormap([0 1 1]);
    figure;plot(sqrt(y(:,1).^2+y(:,2).^2)-R0,y(:,3),'.');xlabel('r');ylabel('Z');
    xlim([-0.5*R0,0.5*R0]);ylim([-0.5*R0,0.5*R0]);axis equal;
elseif(emtype==6)
    title('E=[0,E0,0], B=[0,0,B0], Uniform E&B field');axis equal;
elseif(emtype==7) % non-uniform E
    title('E=[E0*cos(x),0,0], B=[0,0,B0], Non-uniform E');axis equal;
else
    title('E=0, B=[0,0,B0], Uniform B field');axis equal;
end
