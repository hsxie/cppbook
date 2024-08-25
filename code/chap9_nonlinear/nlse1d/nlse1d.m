% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-06-05 12:47
% FD + RK4 to solve NLSE, iE_t + p*E_xx + (V- q*|E|^2)E = 0
% A better method is split-step Fourier method, but only for periodic
% system.
% Test solitons OK.
function nlse1d
    close all; clear; clc;
   
    p=1; q=1;

    L=100; dx=0.001*L; % space grid
    x=(0:dx:L)';
    nx = length(x);
    
    dt = 0.001; % Time steps
    nsteps = floor(10.0/dt);

    % Initial condition
    E = onesoliton(x,0.8,3.2,0.4*L)+onesoliton(x,0.4,-3.2,0.6*L);
% 	E = onesoliton(x,1.8,3.2,0.2*L)+onesoliton(x,2.5,0.5,0.5*L)+...
%         onesoliton(x,1.4,-1.6,0.75*L);
% 	E =1*exp(-(x-0.4*L).^2).*exp(1i*1*x/2);
%     E=sqrt(3)*sech(x-0.4*L).*exp(1i*4*x/2);
    
    figure('units','normalized','position',[0.01 0.12 0.5 0.5],...
        'Color','white','DefaultAxesFontSize',15);
    j=1; t=0;
    for ii=1:nsteps
        
        % Runge-Kutta step
        k1=dt*nlseqn(E,p,q,dx);
        k2=dt*nlseqn(E+k1/2,p,q,dx);
        k3=dt*nlseqn(E+k2/2,p,q,dx);
        k4=dt*nlseqn(E+k3,p,q,dx);
        E=E+k1/6+k2/3+k3/3+k4/6;
        
        % Animate
        if mod(ii,100)==0
            plot(x,real(E),':',x,imag(E),'--',x,abs(E),'LineWidth',2);
            legend('Re(E)','Im(E)','|E|'); legend('boxoff');
            title(['p=',num2str(p),', q=',num2str(q),...
                ', Nx=',num2str(nx),', dt=',num2str(dt),...
                ', t=',num2str(t)]);
            xlabel('x'); ylabel('E');
            axis([0,L,-2,2]);
            drawnow;
            % pause(0.01);
            F(j)=getframe(gcf);
            j=j+1;
            print(gcf,'-dpng',['nlse1d_p=',num2str(p),',q=',num2str(q),...
                ',j=',num2str(j),'.png']);
        end
        t=t+dt;
    end

%     movie(F);
    str=['nlse1d_p=',num2str(p),',q=',num2str(q),'.gif'];
    writegif(str,F,0.1);
    close all;

end

function E=onesoliton(x,a0,v0,x0)
    E=sqrt(2)*a0*sech(a0*(x-x0)).*exp(1i*v0*x/2);
end

function dEdt=nlseqn(E,p,q,dx)
    % NLSE equation: E_t = i*p*E_xx + i*(V- q*|E|^2)E
    V=0;
    E = [E(end-1:end); E; E(1:2)]; % periodic b.c., modify for other b.c.
    dEdt = 1i*p.*(E(4:end-1)-2*E(3:end-2)+E(2:end-3))./dx^2 + ...
        1i.*(V+q.*abs(E(3:end-2)).^2).*E(3:end-2);
end

