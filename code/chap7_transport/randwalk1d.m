% 2013-03-04 09:49
% ref: MIT OpenCourse - Plasma Transport Theory, Fall 2003, ex 1.5
close all; clear; clc;
Np=10000; Nt=1000;

si=[1:1:10,sqrt(24)];
Di_theory=si.^2/24; Di=0.0.*si;
ns=length(si);
for is=1:ns
    xp=zeros(1,Np);
    for it=1:Nt
        Rn=(rand(1,Np)-0.5); % rand number in (-0.5, 0.5)
        dxn=si(is)*Rn;
        xp=xp+dxn;
    end
    xp2_avg=mean(xp.^2);
    Di(is)=xp2_avg/(2*Nt); % diffusion coefficient
end

h=figure('unit','normalized','Position',[0.01 0.17 0.5 0.6]);
set(gcf,'DefaultAxesFontSize',12);

subplot(211);histfit(xp);
[mu,sigma,muci,sigmaci] = normfit(xp);
xlabel('X'); ylabel('P');
title(['s^2=',num2str(si(ns)^2),', N=',num2str(Nt),', <X^2>=',num2str(xp2_avg),...
    ', D=',num2str(Di(ns)),', D_{theory}=',num2str(Di_theory(ns))]);
legend('hist','normfit, P \propto e^{-(x-\mu)^2/\sigma^2}');
legend('boxoff');

subplot(212);plot(si(1:ns-1),Di(1:ns-1),'g-*',si(1:ns-1),Di_theory(1:ns-1),'r--','LineWidth',2);
% hold on; plot(si(ns),Di(ns),'o',si(ns),Di_theory(ns),'*');
title('D=s^2/24'); xlabel('s'); ylabel('D');
legend('Rand Walk','Theory',2); legend('boxoff');

print('-dpng','randwalk1d.png');
