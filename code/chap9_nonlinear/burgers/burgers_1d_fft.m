% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2017-04-05 10:08
% Solve the 1D Burgers equation, spectral method

% close all; clear; clc;
function burgers_1d_fft()
    global nu kx;
    close all;
    nu=0.02; nt=1001; dt=0.02; 
    L=16; N=256*1;
    x=L/N*[-N/2:N/2-1]; dx=x(2)-x(1);
    kx=(2*pi/L)*[0:N/2-1 -N/2:-1].';
    u=exp(-(x+3).^2); ut=fft(u);
    t=(1:nt)*dt;
    [t,utsol]=ode45(@burgers_rhs,t,ut);
    usol=ifft(utsol,[],2);

    figure('unit','normalized','position',[0.1,0.1,0.4,0.5],...
        'DefaultAxesFontSize',12);
    for it=1:floor(nt/5):nt
        if(it<=1)
            plot(x,usol(it,:),'r:','LineWidth',2);hold on;
        else
            plot(x,usol(it,:),'b','LineWidth',2);hold on;
        end
        [ym,idx]=max(usol(it,:));
        text(x(idx),ym+0.05,['t=',num2str(t(it))]);
    end
    title(['Burgers FFT, \nu=',num2str(nu),', L=',num2str(L),', dx=',num2str(dx),', dt=',num2str(dt),...
        ', nt=',num2str(nt)]);
    xlim([min(x),max(x)]); ylim([0,1.1]); xlabel('x');ylabel('u');
    print(gcf,'-dpng',['burgers_fft_nu=',num2str(nu),',L=',num2str(L),...
        ',dx=',num2str(dx),',dt=',num2str(dt),...
        ',nt=',num2str(nt),'.png']);
end

function dut=burgers_rhs(t,ut)
    global nu kx;
    u=ifft(ut);
%     dut=-nu*kx.^2.*ut-1i*kx.*fft(u.*u); % wrong
    dut=-nu*kx.^2.*ut-fft(u.*ifft(1i*kx.*ut));
end
