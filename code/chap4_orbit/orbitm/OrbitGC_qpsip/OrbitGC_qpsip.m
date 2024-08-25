% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2012-12-22 12:22
% OrbitGC.m, Guiding center orbit in Tokamak
% Ref: R. B. White 2001Book p70 and ORBIT code
% epsilon=r/R0<<1, alaytical equilibrium, concentric circles
% time unit: Omega=q*B0/m=1.0, on axis B0=1
% length unit: R0=1.0
% q=q(psip), no E field
% Update: 2012-12-27 23:26
% This version: r=sqrt(2*psi), d(psi)/d(psip)=q(psip)
% psiw=psip at wall (r=a)
% test OK

function OrbitGC_qpsip
    close all; clear; clc;
   
    global q1 q2 q3 mu psiw;
    psiw=0.043636; % psiw=0.043636, a=0.40
    q1=1.0; q2=1.0; q3=1.0;
    
    a=sqrt(2*psiw*(q1+q2/2+q3/3)); % minor radius
    
    % initial values
    psip0=0.6*psiw; theta0=0*pi/4; zeta0=0;
    psin0=psip0/psiw;
    r0=sqrt(2*psip0*(q1+q2/2*psin0+q3/3*psin0^2));
    q=q1+q2*psin0+q3*psin0^2;

    % calculate initial parameters
    R=1+r0*cos(theta0);
    g=1.0; Bt=g/R;
    Bp=r0/(q*R);
    B=sqrt(Bt^2+Bp^2);
    lambda0=0.9; % pitch angle = mu*B/E
    % v ~ rhoi/Omega ~ (10^-3)R0/Omega, then E0 ~ v^2 ~ 10^-6
    E0=1/5e4; % sensitive, E0> e.g., 1/2e6..., banana
    mu=lambda0*E0;
    drc=1; % v0= +/- give different orbits
    v0=sqrt(2*E0);
	rhopara0=drc*v0*sqrt(1-lambda0*B)/B;
    
    % solve
    y0=[zeta0, theta0, psip0, rhopara0];
    options=odeset('RelTol',1e-10,'AbsTol',[1e-10 1e-10 1e-11 1e-10],'MaxStep',0.1);
%     [t,y] = ode45(@orbit,0:200,y0,options);
    tend=150/abs(rhopara0); dt=tend/2e4;
    [t,y] = ode45(@orbit,0:dt:tend,y0);
    
    % plot
    figure; set(gcf,'DefaultAxesFontSize',15);
    zeta=y(:,1); theta=y(:,2); psip=y(:,3); rhopara=y(:,4);
    r=sqrt(2*q.*psip); R=1+r.*cos(theta); x2=r.*cos(theta); y2=r.*sin(theta);
    x3=R.*cos(zeta); y3=R.*sin(zeta); z3=r.*sin(theta);
    subplot(221);plot(t,zeta,'r-','LineWidth',2);
    axis tight; xlabel('t/\Omega_c'); ylabel('\zeta'); 
    title(['\zeta-t',', E=',num2str(E0),', \Lambda=',num2str(lambda0)]);
    subplot(222);plot(t,rhopara,'r-','LineWidth',2); axis tight; 
    xlabel('t/\Omega_c'); ylabel('\rho_{||}'); 
    title(['\rho_{||}-t',', r0=',num2str(r0),',q=',num2str(q)]);
    subplot(223); plot(x2,y2,'LineWidth',2); axis equal;
    xlabel('x'); ylabel('y'); 
    title(['poloidal projection, a=',num2str(a)]); hold on;
    plot(a.*cos(0:pi/20:2*pi),a.*sin(0:pi/20:2*pi),'r--');
    subplot(224); plot3(x3,y3,z3,'g',x3(1),y3(1),z3(1),'r*');    
    title(['direction=',num2str(drc)]);
    axis equal; % title('3D plot'); view(2);
    
    print(gcf,'-dpng',['E=',num2str(E0),',r=',num2str(r0),',Lambda=',...
    num2str(lambda0),',a=',num2str(a),...
    ',drc=',num2str(drc),'.png']);
end

function dy=orbit(t,y)
        
    global q q1 q2 q3 mu psiw;
    
    dy=zeros(4,1);
    zeta=y(1); theta=y(2); psip=y(3); rhopara=y(4);
    psin=psip/psiw; % normalization psip
    q=q1+q2*psin+q3*psin^2;
    
    % field and other parameters
    psi=psip*(q1+q2/2*psin+q3/3*psin^2);
    r=sqrt(2*psi); % r/R0
    R=1+r*cos(theta);
    g=1.0; % g: poloidal current outside psi
    gppsi=0; % p(prime): derivative of psi
    I=r^2/q; % I: toroidal current inside psi
    Ippsi=2/q^3*(q^2-psin*(q2+2*q3*psin)*(q1+q2/2*psin+q3/3*psin^2));
    Bt=g/R;
    Bp=r/(q*R); % Bp=r/(q*R)
    B=sqrt(Bt^2+Bp^2);
    Bppsip=(r*(1-r^2*R*(q2+2*q3*psin)/psiw/q^2)-g^2*q^2*cos(theta))/(B*R^3*r*q); % dB/dpsip
    Bptheta=B*r*sin(theta)/(R); % dB/dtheta
    % Bpzeta=0;
    
    D=g*q+I+rhopara*(g*Ippsi-I*gppsi);

%     E=rhopara^2*B^2/2+mu*B;

    dy(1)=rhopara*B^2*(q+rhopara*Ippsi)/D-(mu+rhopara^2*B)*I*Bppsip/D;
    dy(2)=rhopara*B^2*(1-rhopara*gppsi)/D+(mu+rhopara^2*B)*g*Bppsip/D;
    dy(3)=-g*(mu+rhopara^2*B)*Bptheta/D;
    dy(4)=-(mu+rhopara^2*B)*((1-rhopara*gppsi)*Bptheta)/D;
end