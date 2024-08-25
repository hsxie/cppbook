% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2013-10-17 19:48
% Solve the L-H equation in Kim & Diamond, 2003, PRL
function Kim03
    close all; clear; clc;
    global a1 a2 a3 b1 b2 b3 c1 c2 d;
    a1=0.2; a2=0.7; a3=0.7; b1=1.5; b2=1; b3=1; c1=1; c2=0.5; d=1;
    options = odeset('RelTol',1e-5,'AbsTol',[1e-4 1e-4 1e-5]);
    [t,y] = ode45(@rhs,[0 200],[0.01 0.01 0],options);
    E=y(:,1); Vzf=y(:,2); N=y(:,3); Q=0.01*t;
    
    figure('unit','normalized','Position',[0.01 0.27 0.45 0.5],...
        'DefaultAxesFontSize',15);
    plot(Q,E,'-',Q,Vzf,'-.',Q,N/5,'--','LineWidth',2);
    legend('E','V_{ZF}','N/5');legend('boxoff');
    xlabel('Q=0.01t');
    ylim([0,1.5]);
    str=['a1=',num2str(a1),', a2=',num2str(a2),', a3=',num2str(a3),...
        ', b1=',num2str(b1),', b2=',num2str(b2),', b3=',num2str(b3),...
        ', c1=',num2str(c1),', c2=',num2str(c2),', d=',num2str(d)];
    title(str);
    print('-dpng',['Kim03_',str,'.png']);
    
end

function dy=rhs(t,y)
    dy=zeros(3,1);
    global a1 a2 a3 b1 b2 b3 c1 c2 d;
    E=y(1); Vzf=y(2); N=y(3); V=d*N^2; Q=0.01*t;
    dy(1)=E*N-a1*E^2-a2*V^2*E-a3*Vzf^2*E;
    dy(2)=b1*E*Vzf/(1+b2*V^2)-b3*Vzf;
    dy(3)=-c1*E*N-c2*N+Q;
end
