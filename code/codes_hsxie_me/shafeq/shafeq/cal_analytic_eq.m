% Hua-sheng XIE, IFTS-ZJU, hushengxie@gmail.com, 2012-07-02 16:28
% an auxiliary code for analytical.F90
close all; clear; clc;

% run setpath;


% numeq=2
% psiw=0.02895; % psi on wall
% q1=1.2886;q2=0.8;q3=-0.4;
% ne1=0.473;ne2=0.23;ne3=0.07;
% te1=0.0;te2=0.0;te3=1.0;

% numeq=1, default Cyclone case
psiw=0.0375;
q1=0.82;q2=1.1;q3=1.0;
ne1=0.205;ne2=0.30;ne3=0.4;
te1=0.415;te2=0.18;te3=0.4;

% new Cyclone case
psiw=0.0381;
q1=0.85;q2=1.04;q3=1.0;
ne1=0.206;ne2=0.31;ne3=0.4;
te1=0.4133;te2=0.183;te3=0.4;

% % HL-2A, D. F. Kong's case
% psiw=0.016;
% q1=1.02;q2=1.0;q3=1.0;
% % ne1=0.466;ne2=0.85;ne3=0.04; % eq=11, before 2014-10-25
% % % ne1=0.466;ne2=0.81;ne3=0.10; % eq=11, new0
% % % ne1=0.466;ne2=0.83;ne3=0.08; % eq=11, new1
% % te1=0.466;te2=0.85;te3=0.04;
% 
% % ne1=0.466;ne2=0.81;ne3=0.08; % eq=11, new2
% % te1=0.466;te2=0.81;te3=0.08;
% 
% % ne1=0.466;ne2=0.84;ne3=0.03; % eq=11, new3
% % te1=0.466;te2=0.84;te3=0.03;
% 
% ne1=0.466;ne2=0.82;ne3=0.06; % eq=11, new4
% te1=0.466;te2=0.82;te3=0.06;
% 
% % ne1=0.466;ne2=0.82;ne3=0.05; % eq=11, new5
% % te1=0.466;te2=0.82;te3=0.05;
% 
% 
% % eq=8, xhs ballooning
% % psiw=0.012;
% % q1=1.0; q2=0.0; q3=2.0;
% psiw=0.0134;
% q1=1.0; q2=1.0; q3=0.0;
% % ne1=0.4944; ne2=0.6; ne3=0.1;
% ne1=0.46; ne2=0.5; ne3=0.1;
% te1=0.0; te2=0.0; te3=1.0;

psi=0:0.005*psiw:psiw;
psi_n=psi./psiw; % normalization psi

q=q1+q2.*psi_n+q3.*psi_n.^2;

r=sqrt(2.*(q1.*psi_n+(q2.*psi_n.^2)./2.0+(q3.*psi_n.^3)./3.0)*psiw); % r/R0
a=max(r); % minor radius, a/R0
r_n=r./a; % r/a

ne=1.0+ne1.*(tanh((ne2-(psi_n))./ne3)-1.0);
te=1.0+te1.*(tanh((te2-(psi_n))./te3)-1.0);

r0=165; % major radius, unit=cm
b0=13500; % on-axis magnetic field, unit=gauss
etemp0=1e3; % on-axis electron temperature, unit=ev
eden0=0.15e14; % on-axis electron number density, unit=1/cm^3


R=1+r;
betae0=4.03e-11*etemp0*eden0/(b0*b0);
bb=(1./R).*sqrt(1+(r./(q.*R)).^2); % normalized B, right?
% betae=betae0.*te.*ne./(bb.*bb);   % Q: how to give the B^2 field?
betae=betae0.*te.*ne;   % Q: how to give the B^2 field?

% tau_A/tau_cs0, tau_cs0 is GTC unit, tau_A(r)=q(r)*R/v_A(r) is local Alfven time
% taua_o_taucs0=q.*R.*sqrt(betae./2.0);
taua_o_taucs0=q.*sqrt(betae./2.0);

hf=figure('unit','normalized','Position',[0.1 0.2 0.7 0.7],...
            'Name','Equilibrium profile',... % 'menubar','none',...
            'NumberTitle','off');
set(gcf,'DefaultAxesFontSize',14);

subplot(341);plot(psi_n,q,'--g',r_n,q,'r','LineWidth',2);xlim([0,1]);ylim([0,4]);grid on;
title(['q=/',num2str(q1),',',num2str(q2),',',num2str(q3),'/']);
legend('(psi)','(r)',2);legend('boxoff');

subplot(342);plot(psi_n,ne,'--g',r_n,ne,'r','LineWidth',2);xlim([0,1]);grid minor;
title(['ne=/',num2str(ne1),',',num2str(ne2),',',num2str(ne3),'/']);

subplot(343);plot(psi_n,te,'--g',r_n,te,'r','LineWidth',2);xlim([0,1]);grid on;
title(['Te=/',num2str(te1),',',num2str(te2),',',num2str(te3),'/']);

subplot(344);plot(psi_n,betae,'--g',r_n,betae,'r','LineWidth',2);xlim([0,1]);grid minor;
title(['\beta_e=',num2str(betae0)]);  % betae(1) \neq betae0 !!

subplot(345);
plot(psi_n(1:end-1),(log(q(2:end))-log(q(1:end-1)))./(psi_n(2:end)-psi_n(1:end-1)),...
    '--g',r_n(1:end-1),(log(q(2:end))-log(q(1:end-1)))./(r_n(2:end)-r_n(1:end-1)),'r','LineWidth',2);
xlim([0,1]);title('dln(q)');grid minor;

subplot(346);
dlnne=-(log(ne(2:end))-log(ne(1:end-1)))./(r_n(2:end)-r_n(1:end-1));
plot(psi_n(1:end-1),-(log(ne(2:end))-log(ne(1:end-1)))./(psi_n(2:end)-psi_n(1:end-1)),...
    '--g',r_n(1:end-1),dlnne,'r','LineWidth',2);
xlim([0,1]);title('-dln(ne)');grid on;

subplot(347);
dlnte=-(log(te(2:end))-log(te(1:end-1)))./(r_n(2:end)-r_n(1:end-1));
plot(psi_n(1:end-1),-(log(te(2:end))-log(te(1:end-1)))./(psi_n(2:end)-psi_n(1:end-1)),...
    '--g',r_n(1:end-1),dlnte,'r','LineWidth',2);
xlim([0,1]);title('dln(Te)');grid on;
% plot(psi_n(1:end-1),(te(2:end)-te(1:end-1))./(psi_n(2:end)-psi_n(1:end-1)),...
%     '--g',r_n(1:end-1),(te(2:end)-te(1:end-1))./(r_n(2:end)-r_n(1:end-1)),'r','LineWidth',2);
% xlim([0,1]);title('dTe');grid on;

subplot(348);plot(psi_n,taua_o_taucs0,'--g',r_n,taua_o_taucs0,'r','LineWidth',2);xlim([0,1]);grid minor;
title('\tau_A/\tau_{cs0}');

subplot(349);plot(psi_n,r,'--g','LineWidth',2);
grid on;xlim([0,1]);title(['r(psi), a/R0=',num2str(a)]);

subplot(3,4,10);  % s=(r/q)*(dq/dr)
s=r_n(1:end-1).*(log(q(2:end))-log(q(1:end-1)))./(r_n(2:end)-r_n(1:end-1));
plot(psi_n(1:end-1),s,'--g',r_n(1:end-1),s,'r','LineWidth',2);
xlim([0,1]);title('shear');grid minor;

subplot(3,4,11);  % alpha=-q^2*R*(dbeta/dr), Q: how to give R?
alpha=-q(1:end-1).^2.*(R(1:end-1).*2).*(betae(2:end)-betae(1:end-1))./(r(2:end)-r(1:end-1));
plot(psi_n(1:end-1),alpha,'--g',r_n(1:end-1),alpha,'r','LineWidth',2);
xlim([0,1]);title('\alpha');grid minor;


print(hf,'-dpng','analytical_equilibrium_profile.png');

%% 
riflux=0.5;
a
psiiflux=interp1(r_n,psi_n,riflux,'spline')
qiflux=interp1(r_n,q,riflux,'spline')
siflux=interp1(r_n(1:end-1),s,riflux,'spline')
teiflux=interp1(r_n,te,riflux,'spline')/te(1)
neiflux=interp1(r_n,ne,riflux,'spline')/ne(1)
Lneiflux=interp1(r_n(1:end-1),dlnne,riflux,'spline')/a
Lteiflux=interp1(r_n(1:end-1),dlnte,riflux,'spline')/a

r0iflux=interp1(psi_n,r_n,0.3,'spline')
r1iflux=interp1(psi_n,r_n,0.92,'spline')
r0iflux+(r1iflux-r0iflux)*0.6

%%
% dat=load([path,'gtc_rg.out']);
% 
% tmp=interp1(r_n,psi_n,dat(:,2),'spline')
% subplot(221);plot(dat(:,2));
% subplot(222);plot(diff(dat(:,2)));
% subplot(223);plot(tmp);
% subplot(224);plot(diff(tmp));
