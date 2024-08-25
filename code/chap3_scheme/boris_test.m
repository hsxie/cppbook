% Hua-sheng XIE, 2015-05-10 23:58
close all; clear; clc;
E0=[0,-0.4,0]; B0=[0,0,1.5]; tmp=1; nt=1000*tmp; dt=0.02/tmp; qm=1.0;
Ex=E0(1); Ey=E0(2); Ez=E0(3); Bx=B0(1); By=B0(2); Bz=B0(3);
vx(1)=0.8; vy(1)=0.1; vz(1)=0.2; x(1)=0.0; y(1)=0.0; z(1)=0.0;
vperp=sqrt(vx(1)^2+vy(1)^2); vd=[Ey/Bz,0,vz(1)]; wc=qm*Bz; rL=vperp/(qm*Bz); 
rc=[rL*vy(1)/vperp+x(1),-rL*vx(1)/vperp+y(1),z(1)]; phi0=atan(-vy(1)/vx(1));
t=(1:nt)*dt; rcx=rc(1)+vd(1)*t; rcy=rc(2)+vd(2)*t; rcz=rc(3)+vd(3)*t;
% need add the correct form of exact solution
% rx=rcx+rL*sin(pi/2)*sin(wc*t+phi0); ry=rcy+rL*cos(wc*t+phi0); rz=rcz;
figure('unit','normalized','position',[0.02,0.1,0.6,0.4],'DefaultAxesFontSize',15);
for method=1:2
    for it=1:nt
        if(method==2)
            vx(it+1)=vx(it)+qm*(Ex+(vy(it)*Bz-vz(it)*By))*dt;
            vy(it+1)=vy(it)+qm*(Ey+(vz(it)*Bx-vx(it)*Bz))*dt;
            vz(it+1)=vz(it)+qm*(Ez+(vx(it)*By-vy(it)*Bx))*dt;
        else
            qtmp=dt*qm/2;
            hx=qtmp*Bx; hy=qtmp*By; hz=qtmp*Bz; h2=hx*hx+hy*hy+hz*hz;
            sx=2*hx/(1+h2); sy=2*hy/(1+h2); sz=2*hz/(1+h2);
            ux=vx(it)+qtmp*Ex; uy=vy(it)+qtmp*Ey; uz=vz(it)+qtmp*Ez;
            uxtmp=ux+(uy*sz-uz*sy)+((uz*hx-ux*hz)*sz-(ux*hy-uy*hx)*sy);
            uytmp=uy+(uz*sx-ux*sz)+((ux*hy-uy*hx)*sx-(uy*hz-uz*hy)*sz);
            uztmp=uz+(ux*sy-uy*sx)+((uy*hz-uz*hy)*sy-(uz*hx-ux*hz)*sx);
            vx(it+1)=uxtmp+qtmp*Ex;
            vy(it+1)=uytmp+qtmp*Ey;
            vz(it+1)=uztmp+qtmp*Ez;
        end
        x(it+1)=x(it)+vx(it+1)*dt;
        y(it+1)=y(it)+vy(it+1)*dt;
        z(it+1)=z(it)+vz(it+1)*dt;
    end    
    if(method==1)
        subplot(121); plot3(x,y,z,'b-','Linewidth',2); hold on;
        subplot(122); plot3(vx,vy,vz,'b-','Linewidth',2); hold on;
    else
        subplot(121); plot3(x,y,z,'g--','Linewidth',2); hold on;
        subplot(122); plot3(vx,vy,vz,'g--','Linewidth',2); hold on;
    end
end
subplot(121); box on; xlabel('x'); ylabel('y'); zlabel('z');
hold on; plot3(rcx,rcy,rcz,'r--','Linewidth',2); axis equal;
% hold on; plot3(rx,ry,rz,'m--','Linewidth',2);
legend('Boris','Euler','Guiding center',1); legend('boxoff');
subplot(122); box on; xlabel('v_x'); ylabel('v_y'); zlabel('v_z');
