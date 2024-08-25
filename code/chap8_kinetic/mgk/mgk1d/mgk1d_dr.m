% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2017-02-04 15:39
% Solve the ITG integral equation in ballooning space, Ref: Dong1992
% Ref also Zhi-xin LU's version
% 17-02-16 12:52
close all; clear; clc;

global tau kt q s epsn etai tk;
global dvx dvy vxx vyy dth tt;
global nL hll Mij;

runtime=cputime;

dat=[0.45    -0.615095+0.263450i;
    0.50    -0.619550+0.245009i;
    0.55    -0.623357+0.222611i]; % From HanMK, Dong92
id=1;
kt=dat(id,1)/sqrt(2); % k_theta*rho_i
wr=real(dat(id,2)); wi=imag(dat(id,2));
s=1.0; q=1.0; tau=1.0; epsn=0.25; etai=2.5;
tk=0.0; % theta_k

wt0=(wr+1i*wi);
wt=wt0*kt/epsn;

nth=400; thmax=pi*4.0; thmin=-thmax;
dth=(thmax-thmin)/(nth-1); tt=thmin:dth:thmax; %nt=length(tt);

nvx=25; nvy=25; vxmax=5.0; vymax=5.0; dvx=vxmax/nvx; dvy=vymax/nvy;
vxx=(-0.5+(1:nvx))*dvx; vyy=(-0.5+(1:nvy))*dvy; % be careful to avoid vpara=0

nL=30; % # of Hermite basis functions
wg=wt; % initial guess

c1=0.7; tc=0.0; tt1=(tt-tc)*c1;
hll=zeros(nL,length(tt)); 
for l=0:(nL-1)
    hl=exp(-tt1.^2/2).*hermiteH(l,tt1)*sqrt(c1/(sqrt(pi)*(2^l)*factorial(l)));
	hll(l+1,:)=hl;
end

detM=@(w)fun_detM(w);
[x, fval]=fsolve(detM,wg);

[v,d]=eig(Mij);
phi=v(:,1).'*hll;
phi=phi./phi(floor(nth/2));
% phi=phi./phi(floor(nth/2)+1);

runtime=cputime-runtime;
%%
close all;
h = figure('Unit','Normalized','position',...
    [0.02 0.13 0.4 0.5],'DefaultAxesFontSize',15);
plot(tt,real(phi),'b',tt,imag(phi),'r--','linewidth',2);
legend('Re','Im','abs');ylabel('\delta\phi'); legend('boxoff');

th_dr=tt;
phi_dr=phi;

if(id==1)
run phi_hd7_dat; hold on;
plot(th_hd7,real(phi_hd7),'k:',th_hd7,imag(phi_hd7),'m:','Linewidth',2);
end

title(['\omega=',num2str(x,4),', runtime=',num2str(runtime),'s']);
xlabel(['\theta, \omega^{T}=',num2str(wt,4)]);
print(gcf,'-dpng',['mgk1d_dr_s=',num2str(s),',q=',num2str(q),...
    ',tau=',num2str(tau),',eta=',num2str(etai),',epsn=',num2str(epsn),...
    ',kt=',num2str(kt),',nL=',num2str(nL),',nth=',num2str(nth),...
    ',nvx=',num2str(nvx),',nvy=',num2str(nvy),...
    ',thmax=',num2str(thmax),',vxmax=',num2str(vxmax),...
    ',vymax=',num2str(vymax),'.png']);

