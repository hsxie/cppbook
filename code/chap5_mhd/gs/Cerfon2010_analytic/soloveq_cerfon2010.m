% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-01-08 09:39
% Equilibrium Solver, Analytic Solution, Solovev type
% Ref: [Cerfon2010] Cerfon, A. J. & Freidberg, J. P., ¡°One size fits all¡±
%      analytic solutions to the Grad¨CShafranov equation, Phys. Plasmas,
%      2010, 17, 032502.
% Elongated "D" shape: x=1+epsilon*cos(tau+alpha*sin(tau)), R=R0*x
%                      y=epsilon*kappa*sin(tau), Z=R0*y
% epsilon=a/R0, kappa is elongation, sin(alpha)=delta is triangularity
% 
% inull: =0, no null;  =1 single null; =2 double null
% Test OK, 2013-01-11 01:18
function soloveq_cerfon2010
close all; clear; clc;
    global fBR fBZ fBp fBt;
inull=0;
epsilon=0.32;
kappa=1.7;
delta=0.33;
alpha=asin(delta);
N1=-(1+alpha)^2/(epsilon*kappa^2);N2=(1-alpha)^2/(epsilon*kappa^2);
N3=-kappa/(epsilon*cos(alpha)^2);
A=-0.155;
xsep=1-1.1*delta*epsilon;
ysep=-1.1*kappa*epsilon;

Psi0=1.0; % normalization psi
B0=1.0; % on axis B field

%% psi, 0-7 (even symmetrtric), 8-12 (odd symmetrtric, for single null)
psi0=@(x,y)x.^4/8+A.*(x.*log(x)/2-x.^4/8);
% a1=1;a2=-1;psi0=@(x,y)a1*x.^4/8-a2*y.^2/2; % SunYW's psi_P
psi1=@(x,y)1;
psi2=@(x,y)x.^2;
psi3=@(x,y)y.^2-x.^2.*log(x);
psi4=@(x,y)x.^4-4*x.^2.*y.^2;
psi5=@(x,y)2*y.^4-9*y.^2.*x.^2+3*x.^4.*log(x)-12*x.^2.*y.^2.*log(x);
psi6=@(x,y)x.^6-12*x.^4.*y.^2+8*x.^2.*y.^4;
psi7=@(x,y)8*y.^6-140*y.^4.*x.^2+75*y.^2.*x.^4-15*x.^6.*log(x)+...
    180*x.^4.*y.^2.*log(x)-120*x.^2.*y.^4.*log(x);
psi8=@(x,y)y;
psi9=@(x,y)y.*x.^2;
psi10=@(x,y)y.^3-3*y.*x.^2.*log(x);
psi11=@(x,y)3*y.*x.^4-4*x.^2.*y.^3;
% psi11=@(x,y)3*y.*x.^4-4*x.^2.*y.^5; % SunYW's psi11
psi12=@(x,y)8*y.^5-45*y.*x.^4-80*x.^2.*y.^3.*log(x)+60*y.*x.^4.*log(x);

%% psix, psixx, psiy, psiyy, 0-7, 8-12
psi0x=@(x,y)A.*(log(x)./2-x.^3./2+1./2)+x.^3./2;
psi0xx=@(x,y)A.*(1./(2.*x)-(3.*x.^2)./2)+(3.*x.^2)./2;
psi0y=@(x,y)0;
psi0yy=@(x,y)0;
% psi0x=@(x,y)a1.*x.^3/2; % SunYW's psi_P
% psi0xx=@(x,y)a1.*(3.*x.^2)./2;
% psi0y=@(x,y)-a2.*y;
% psi0yy=@(x,y)-a2;

psi1x=@(x,y)0;
psi1xx=@(x,y)0;
psi1y=@(x,y)0;
psi1yy=@(x,y)0;

psi2x=@(x,y)2.*x;
psi2xx=@(x,y)2;
psi2y=@(x,y)0;
psi2yy=@(x,y)0;

psi3x=@(x,y)-x-2.*x.*log(x);
psi3xx=@(x,y)-2.*log(x)-3;
psi3y=@(x,y)2.*y;
psi3yy=@(x,y)2;

psi4x=@(x,y)4.*x.^3-8.*x.*y.^2;
psi4xx=@(x,y)12.*x.^2-8.*y.^2;
psi4y=@(x,y)-8.*x.^2.*y;
psi4yy=@(x,y)-8.*x.^2;

psi5x=@(x,y)12.*x.^3.*log(x)-30.*x.*y.^2+3.*x.^3-24.*x.*y.^2.*log(x);
psi5xx=@(x,y)36.*x.^2.*log(x)-24.*y.^2.*log(x)+21.*x.^2-54.*y.^2;
psi5y=@(x,y)8.*y.^3-18.*x.^2.*y-24.*x.^2.*y.*log(x);
psi5yy=@(x,y)24.*y.^2-18.*x.^2-24.*x.^2.*log(x);

psi6x=@(x,y)6.*x.^5-48.*x.^3.*y.^2+16.*x.*y.^4;
psi6xx=@(x,y)30.*x.^4-144.*x.^2.*y.^2+16.*y.^4;
psi6y=@(x,y)-24.*x.^4.*y+32.*x.^2.*y.^3;
psi6yy=@(x,y)-24.*x.^4+96.*x.^2.*y.^2;

psi7x=@(x,y)480.*x.^3.*y.^2-90.*x.^5.*log(x)-400.*x.*y.^4-15.*x.^5-240.*x.*y.^4.*log(x)+720.*x.^3.*y.^2.*log(x);
psi7xx=@(x,y)2160.*x.^2.*y.^2-450.*x.^4.*log(x)-240.*y.^4.*log(x)-165.*x.^4-640.*y.^4+2160.*x.^2.*y.^2.*log(x);
psi7y=@(x,y)150.*x.^4.*y-560.*x.^2.*y.^3+48.*y.^5+360.*x.^4.*y.*log(x)-480.*x.^2.*y.^3.*log(x);
psi7yy=@(x,y)360.*x.^4.*log(x)-1680.*x.^2.*y.^2+150.*x.^4+240.*y.^4-1440.*x.^2.*y.^2.*log(x);

psi8x=@(x,y)0;
psi8xx=@(x,y)0;
psi8y=@(x,y)1;
psi8yy=@(x,y)0;

psi9x=@(x,y)2.*x.*y;
psi9xx=@(x,y)2.*y;
psi9y=@(x,y)x.^2;
psi9yy=@(x,y)0;

psi10x=@(x,y)-3.*x.*y-6.*x.*y.*log(x);
psi10xx=@(x,y)-9.*y-6.*y.*log(x);
psi10y=@(x,y)3.*y.^2-3.*x.^2.*log(x);
psi10yy=@(x,y)6.*y;

psi11x=@(x,y)12.*x.^3.*y-8.*x.*y.^3;
psi11xx=@(x,y)36.*x.^2.*y-8.*y.^3;
psi11y=@(x,y)3.*x.^4-12.*x.^2.*y.^2;
psi11yy=@(x,y)-24.*x.^2.*y;
% psi11x=@(x,y)12.*x.^3..*y-8.*x.*y.^5; % SunYW's psi11
% psi11xx=@(x,y)36.*x.^2..*y-8.*y.^5;
% psi11y=@(x,y)3.*x.^4-20.*x.^2.*y.^4;
% psi11yy=@(x,y)-80.*x.^2.*y.^3;

psi12x=@(x,y)240.*x.^3.*y.*log(x)-120.*x.^3.*y-160.*x.*y.^3.*log(x)-80.*x.*y.^3;
psi12xx=@(x,y)720.*x.^2.*y.*log(x)-120.*x.^2.*y-240.*y.^3-160.*y.^3.*log(x);
psi12y=@(x,y)60.*x.^4.*log(x)-45.*x.^4+40.*y.^4-240.*x.^2.*y.^2.*log(x);
psi12yy=@(x,y)160.*y.^3-480.*x.^2.*y.*log(x);


%% for b.c. constraints
c5=0;c6=0;c7=0;c8=0;c9=0;c10=0;c11=0;c12=0;
if(inull==1) % single null, inull==1
    x1=1+epsilon; y1=0; % outer equatorial point
    x2=1-epsilon; y2=0; % inner equatorial point
    x3=1-delta*epsilon; y3=kappa*epsilon; % upper high point
    x4=xsep; y4=ysep; % lower X point
    x5=1+epsilon; y5=0; % outer equatorial point up-down symmetry
    x6=1-epsilon; y6=0; % inner equatorial point up-down symmetry
    x7=1-delta*epsilon; y7=kappa*epsilon; % upper high point maximum
    x8=xsep; y8=ysep; % By=0 at lower X point
    x9=xsep; y9=ysep; % Bx=0 at lower X point
    x10=1+epsilon; y10=0; % outer equatorial point curvature
    x11=1-epsilon; y11=0; % inner equatorial point curvature
    x12=1-delta*epsilon; y12=kappa*epsilon; % high point curvature
    
    Mc=[psi1(x1,y1),psi2(x1,y1),psi3(x1,y1),psi4(x1,y1),psi5(x1,y1),...
        psi6(x1,y1),psi7(x1,y1),psi8(x1,y1),psi9(x1,y1),psi10(x1,y1),...
        psi11(x1,y1),psi12(x1,y1);
        psi1(x2,y2),psi2(x2,y2),psi3(x2,y2),psi4(x2,y2),psi5(x2,y2),...
        psi6(x2,y2),psi7(x2,y2),psi8(x2,y2),psi9(x2,y2),psi10(x2,y2),...
        psi11(x2,y2),psi12(x2,y2);
        psi1(x3,y3),psi2(x3,y3),psi3(x3,y3),psi4(x3,y3),psi5(x3,y3),...
        psi6(x3,y3),psi7(x3,y3),psi8(x3,y3),psi9(x3,y3),psi10(x3,y3),...
        psi11(x3,y3),psi12(x3,y3);
        psi1(x4,y4),psi2(x4,y4),psi3(x4,y4),psi4(x4,y4),psi5(x4,y4),...
        psi6(x4,y4),psi7(x4,y4),psi8(x4,y4),psi9(x4,y4),psi10(x4,y4),...
        psi11(x4,y4),psi12(x4,y4);
        psi1y(x5,y5),psi2y(x5,y5),psi3y(x5,y5),psi4y(x5,y5),...
        psi5y(x5,y5),psi6y(x5,y5),psi7y(x5,y5),psi8y(x5,y5),...
        psi9y(x5,y5),psi10y(x5,y5),psi11y(x5,y5),psi12y(x5,y5);
        psi1y(x6,y6),psi2y(x6,y6),psi3y(x6,y6),psi4y(x6,y6),...
        psi5y(x6,y6),psi6y(x6,y6),psi7y(x6,y6),psi8y(x6,y6),...
        psi9y(x6,y6),psi10y(x6,y6),psi11y(x6,y6),psi12y(x6,y6);
        psi1x(x7,y7),psi2x(x7,y7),psi3x(x7,y7),psi4x(x7,y7),...
        psi5x(x7,y7),psi6x(x7,y7),psi7x(x7,y7),psi8x(x7,y7),...
        psi9x(x7,y7),psi10x(x7,y7),psi11x(x7,y7),psi12x(x7,y7);
        psi1x(x8,y8),psi2x(x8,y8),psi3x(x8,y8),psi4x(x8,y8),...
        psi5x(x8,y8),psi6x(x8,y8),psi7x(x8,y8),psi8x(x8,y8),...
        psi9x(x8,y8),psi10x(x8,y8),psi11x(x8,y8),psi12x(x8,y8);
        psi1y(x9,y9),psi2y(x9,y9),psi3y(x9,y9),psi4y(x9,y9),...
        psi5y(x9,y9),psi6y(x9,y9),psi7y(x9,y9),psi8y(x9,y9),...
        psi9y(x9,y9),psi10y(x9,y9),psi11y(x9,y9),psi12y(x9,y9);
 psi1yy(x10,y10)+N1*psi1x(x10,y10),psi2yy(x10,y10)+N1*psi2x(x10,y10),...
 psi3yy(x10,y10)+N1*psi3x(x10,y10),psi4yy(x10,y10)+N1*psi4x(x10,y10),...
 psi5yy(x10,y10)+N1*psi5x(x10,y10),psi6yy(x10,y10)+N1*psi6x(x10,y10),...
 psi7yy(x10,y10)+N1*psi7x(x10,y10),psi8yy(x10,y10)+N1*psi8x(x10,y10),...
 psi9yy(x10,y10)+N1*psi9x(x10,y10),psi10yy(x10,y10)+N1*psi10x(x10,y10),...
 psi11yy(x10,y10)+N1*psi11x(x10,y10),psi12yy(x10,y10)+N1*psi12x(x10,y10);
 psi1yy(x11,y11)+N2*psi1x(x11,y11),psi2yy(x11,y11)+N2*psi2x(x11,y11),...
 psi3yy(x11,y11)+N2*psi3x(x11,y11),psi4yy(x11,y11)+N2*psi4x(x11,y11),...
 psi5yy(x11,y11)+N2*psi5x(x11,y11),psi6yy(x11,y11)+N2*psi6x(x11,y11),...
 psi7yy(x11,y11)+N2*psi7x(x11,y11),psi8yy(x11,y11)+N2*psi8x(x11,y11),...
 psi9yy(x11,y11)+N2*psi9x(x11,y11),psi10yy(x11,y11)+N2*psi10x(x11,y11),...
 psi11yy(x11,y11)+N2*psi11x(x11,y11),psi12yy(x11,y11)+N2*psi12x(x11,y11);
 psi1xx(x12,y12)+N3*psi1y(x12,y12),psi2xx(x12,y12)+N3*psi2y(x12,y12),...
 psi3xx(x12,y12)+N3*psi3y(x12,y12),psi4xx(x12,y12)+N3*psi4y(x12,y12),...
 psi5xx(x12,y12)+N3*psi5y(x12,y12),psi6xx(x12,y12)+N3*psi6y(x12,y12),...
 psi7xx(x12,y12)+N3*psi7y(x12,y12),psi8xx(x12,y12)+N3*psi8y(x12,y12),...
 psi9xx(x12,y12)+N3*psi9y(x12,y12),psi10xx(x12,y12)+N3*psi10y(x12,y12),...
 psi11xx(x12,y12)+N3*psi11y(x12,y12),psi12xx(x12,y12)+N3*psi12y(x12,y12)
        ];
    Bc=-[psi0(x1,y1); % psi(x1,y1)=0
        psi0(x2,y2); % psi(x2,y2)=0
        psi0(x3,y3); % psi(x3,y3)=0
        psi0(x4,y4); % psi(x4,y4)=0
        psi0y(x5,y5); % psiy(x5,y5)=0
        psi0y(x6,y6); % psiy(x6,y6)=0
        psi0x(x7,y7); % psix(x7,y7)=0
        psi0x(x8,y8); % psix(x8,y8)=0
        psi0y(x9,y9); % psiy(x9,y9)=0
    psi0yy(x10,y10)+N1*psi0x(x10,y10); % psiyy(x10,y10)=-N1*psix(x10,y10)
    psi0yy(x11,y11)+N2*psi0x(x11,y11); % psiyy(x11,y11)=-N2*psix(x11,y11)
    psi0xx(x12,y12)+N3*psi0y(x12,y12)]; % psixx(x12,y12)=-N3*psiy(x12,y12)

    C=Mc\Bc;
    c1=C(1);c2=C(2);c3=C(3);c4=C(4);c5=C(5);c6=C(6);c7=C(7);
    c8=C(8);c9=C(9);c10=C(10);c11=C(11);c12=C(12);
elseif(inull==2) % double null, inull=2    
    x1=1+epsilon; y1=0; % outer equatorial point
    x2=1-epsilon; y2=0; % inner equatorial point
    x3=xsep; y3=ysep; % high point
    x4=xsep; y4=ysep; % By=0 at high point
    x5=xsep; y5=ysep; % Bx=0 at high point
    x6=1+epsilon; y6=0; % outer equatorial point curvature
    x7=1-epsilon; y7=0; % inner equatorial point curvature

    Mc=[psi1(x1,y1),psi2(x1,y1),psi3(x1,y1),psi4(x1,y1),psi5(x1,y1),...
        psi6(x1,y1),psi7(x1,y1);
        psi1(x2,y2),psi2(x2,y2),psi3(x2,y2),psi4(x2,y2),psi5(x2,y2),...
        psi6(x2,y2),psi7(x2,y2);
        psi1(x3,y3),psi2(x3,y3),psi3(x3,y3),psi4(x3,y3),psi5(x3,y3),...
        psi6(x3,y3),psi7(x3,y3);
        psi1x(x4,y4),psi2x(x4,y4),psi3x(x4,y4),psi4x(x4,y4),...
        psi5x(x4,y4),psi6x(x4,y4),psi7x(x4,y4);
        psi1y(x5,y5),psi2y(x5,y5),psi3y(x5,y5),psi4y(x5,y5),...
        psi5y(x5,y5),psi6y(x5,y5),psi7y(x5,y5);
        psi1yy(x6,y6)+N1*psi1x(x6,y6),psi2yy(x6,y6)+N1*psi2x(x6,y6),...
        psi3yy(x6,y6)+N1*psi3x(x6,y6),psi4yy(x6,y6)+N1*psi4x(x6,y6),...
        psi5yy(x6,y6)+N1*psi5x(x6,y6),psi6yy(x6,y6)+N1*psi6x(x6,y6),...
        psi7yy(x6,y6)+N1*psi7x(x6,y6);
        psi1yy(x7,y7)+N2*psi1x(x7,y7),psi2yy(x7,y7)+N2*psi2x(x7,y7),...
        psi3yy(x7,y7)+N2*psi3x(x7,y7),psi4yy(x7,y7)+N2*psi4x(x7,y7),...
        psi5yy(x7,y7)+N2*psi5x(x7,y7),psi6yy(x7,y7)+N2*psi6x(x7,y7),...
        psi7yy(x7,y7)+N2*psi7x(x7,y7);
        ];
    Bc=-[psi0(x1,y1); % psi(x1,y1)=0
        psi0(x2,y2); % psi(x2,y2)=0
        psi0(x3,y3); % psi(x3,y3)=0
        psi0x(x4,y4); % psix(x4,y4)=0
        psi0y(x5,y5); % psiy(x5,y5)=0
        psi0yy(x6,y6)+N1*psi0x(x6,y6); % psiyy(x6,y6)=-N1*psix(x6,y6)
        psi0yy(x7,y7)+N2*psi0x(x7,y7)]; % psiyy(x7,y7)=-N2*psix(x7,y7)
    C=Mc\Bc;
    c1=C(1);c2=C(2);c3=C(3);c4=C(4);c5=C(5);c6=C(6);c7=C(7);
else % no null, inull=0
    x1=1+epsilon; y1=0; % outer equatorial point
    x2=1-epsilon; y2=0; % inner equatorial point
    x3=1-delta*epsilon; y3=kappa*epsilon; % high point
    x4=1-delta*epsilon; y4=kappa*epsilon; % high point maximum
    
    % Mc*C=Bc
    Mc=[psi1(x1,y1),psi2(x1,y1),psi3(x1,y1),psi4(x1,y1);
        psi1(x2,y2),psi2(x2,y2),psi3(x2,y2),psi4(x2,y2);
        psi1(x3,y3),psi2(x3,y3),psi3(x3,y3),psi4(x3,y3);
        psi1x(x4,y4),psi2x(x4,y4),psi3x(x4,y4),psi4x(x4,y4)];
    Bc=-[psi0(x1,y1); % psi(x1,y1)=0
        psi0(x2,y2); % psi(x2,y2)=0
        psi0(x3,y3); % psi(x3,y3)=0
        psi0x(x4,y4)]; % psix(x4,y4)=0
    
    C=Mc\Bc; % matlab: x=A\B for A*x=B; x=B/A for x*A=B
    c1=C(1);c2=C(2);c3=C(3);c4=C(4);
end

%% the Solovev solution psi(x,y)
fPsi=@(x,y)psi0(x,y)+c1.*psi1(x,y)+c2*psi2(x,y)+...
  c3*psi3(x,y)+c4*psi4(x,y)+c5*psi5(x,y)+c6*psi6(x,y)+c7*psi7(x,y)+...
  c8*psi8(x,y)+c9*psi9(x,y)+c10*psi10(x,y)+c11*psi11(x,y)+c12*psi12(x,y);
fPsix=@(x,y)psi0x(x,y)+c1.*psi1x(x,y)+c2*psi2x(x,y)+...
  c3*psi3x(x,y)+c4*psi4x(x,y)+c5*psi5x(x,y)+c6*psi6x(x,y)+...
  c7*psi7x(x,y)+c8*psi8x(x,y)+c9*psi9x(x,y)+c10*psi10x(x,y)+...
  c11*psi11x(x,y)+c12*psi12x(x,y);
fPsiy=@(x,y)psi0y(x,y)+c1.*psi1y(x,y)+c2*psi2y(x,y)+...
  c3*psi3y(x,y)+c4*psi4y(x,y)+c5*psi5y(x,y)+c6*psi6y(x,y)+...
  c7*psi7y(x,y)+c8*psi8y(x,y)+c9*psi9y(x,y)+c10*psi10y(x,y)+...
  c11*psi11y(x,y)+c12*psi12y(x,y);
fPsixx=@(x,y)psi0xx(x,y)+c1.*psi1xx(x,y)+c2*psi2xx(x,y)+...
  c3*psi3xx(x,y)+c4*psi4xx(x,y)+c5*psi5xx(x,y)+c6*psi6xx(x,y)+...
  c7*psi7xx(x,y)+c8*psi8xx(x,y)+c9*psi9xx(x,y)+c10*psi10xx(x,y)+...
  c11*psi11xx(x,y)+c12*psi12xx(x,y);
fPsiyy=@(x,y)psi0yy(x,y)+c1.*psi1yy(x,y)+c2*psi2yy(x,y)+...
  c3*psi3yy(x,y)+c4*psi4yy(x,y)+c5*psi5yy(x,y)+c6*psi6yy(x,y)+...
  c7*psi7yy(x,y)+c8*psi8yy(x,y)+c9*psi9yy(x,y)+c10*psi10yy(x,y)+...
  c11*psi11yy(x,y)+c12*psi12yy(x,y);

% cal Pressure, Bt, BR, BZ, Bp profile
fP=@(x,y)-Psi0^2*(1-A).*fPsi(x,y);
fBt=@(x,y)(B0^2-2*Psi0^2*A.*fPsi(x,y))./x;
fBR=@(x,y)-fPsiy(x,y)./x; % B_R=-(1/R)dPsi/dZ
fBZ=@(x,y)fPsix(x,y)./x; % B_Z=(1/R)dPsi/dR
fBp=@(x,y)sqrt(fBR(x,y).^2+fBZ(x,y).^2);
fB=@(x,y)sqrt(fBR(x,y).^2+fBZ(x,y).^2+fBt(x,y).^2);


%% plot
% figure; set(gcf,'DefaultAxesFontSize',15);
% xmin=0.8*(1-epsilon); xmax=1.2*(1+epsilon); dx=0.01;
% ymin=-1.2*kappa*epsilon; ymax=1.2*kappa*epsilon; dy=0.01;
% [X,Y]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
% Psi=fpsi(X,Y);
% psi_mag=interp2(X,Y,Psi,1,0,'cubic'); psi_b=0;
% cc=linspace(psi_mag,psi_b+0.5*(psi_b-psi_mag),20);
% contour(X,Y,Psi,cc);axis equal;
% title(['inull=',num2str(inull)]);
% tau=0:pi/50:2*pi;
% xx=1+epsilon.*cos(tau+alpha.*sin(tau));
% yy=kappa*epsilon.*sin(tau);
% hold on; plot(xx,yy,'r--','LineWidth',2);
% 
% print('-dpng',['inull=',num2str(inull),'.png']);

%% plot, copy from gseq_Helander2001.m
R0=1; Rx=0.6*R0;

dR=0.005*R0;  % Field line sensitive?
dZ=dR; Rmin=0.5*R0; Rmax=R0+Rmin; Zmax=0.9*R0;Zmin=-Zmax;
[R,Z]=meshgrid(Rmin:dR:Rmax,Zmin:dZ:Zmax);

Psi=fPsi(R,Z); BR=fBR(R,Z); BZ=fBZ(R,Z); Bt=fBt(R,Z); B=fB(R,Z);
Bp=sqrt(BR.^2+BZ.^2); r=sqrt((abs(R)-R0).^2+Z.^2);
Psitmp=fPsi(Rx+0.03*R0,0); % 
ind=find(Psi>Psitmp);
Psiplt=Psi; Psiplt((ind))=NaN; % filter outer separatrix points
Bpplt=Bp; Bpplt((ind))=NaN;
Btplt=Bt; Btplt((ind))=NaN;
Bplt=B; Bplt((ind))=NaN;

indo=find(Psi==min(min(Psi))); % O point, i.e., the axis
Psio=Psi(indo); Ro=R(indo); Zo=Z(indo);
indxp=find(Bp<=eps & r>5*dR); % X point
Psixp=Psi(indxp); Rxp=R(indxp); Zxp=Z(indxp);

subplot(221);contour(R,Z,Psiplt);axis equal; 
title(['\psi(R,Z), Rx=',num2str(Rx/R0),'R0']);
hold on; plot(Ro,Zo,'r+',Rxp,Zxp,'rx');

subplot(222);pcolor(R,Z,Psiplt);shading interp; axis equal; colorbar; 
% title(['\psi(R,Z), \kappa=',num2str(kappa),', \tau=',num2str(tau)]);
hold on; plot(Ro,Zo,'g+',Rxp,Zxp,'gx');

 % for 3D and field line plot
subplot(223);
Ra0=Rx+0.05*R0; Za0=0; % change Ra0
ya0=[Ra0,0,Za0];
[lpa,Ya]=ode45(@fieldline,0:dR:10*R0,ya0); % analytic field line
Rla=Ya(:,1);phila=Ya(:,2);Zla=Ya(:,3);
xla=Rla.*cos(phila);yla=Rla.*sin(phila);zla=Zla;
nlp=length(xla); ipl=floor(nlp/1); % adjust plot range
plot3(xla(1:ipl),yla(1:ipl),zla(1:ipl),'y','LineWidth',2);
title(['Field Line, Ra0=',num2str(Ra0),'R0, B_0=',num2str(B0)]);

extrMaxIndex=find(diff(sign(diff(Zla)))==-2)+1;
ind1=extrMaxIndex(1);ind2=extrMaxIndex(2)+1;
Ra=Rla(ind1:ind2); Za=Zla(ind1:ind2);
Phia=-0.5*pi:pi/20:2*pi/2;
[Rl,Phil]=meshgrid(Ra,Phia);
[Zl,Phil]=meshgrid(Za,Phia);
xl=Rl.*cos(Phil);yl=Rl.*sin(Phil); zl=Zl;
Fl=fBt(Rl,Zl);
hold on; surf(xl,yl,zl,Fl,'EdgeColor','none');axis equal; % alpha(0.5);

% cal and plot q(psi)
nqp=20;
for iq=1:nqp
    Zq=0;Rq=R0-0.8*(iq/nqp)*(R0-Rx);phiq=0;
    yq0=[Rq,phiq,Zq];
    [lpq,Yq]=ode45(@fieldline,0:0.001*R0:10*R0,yq0);
    Rlq=Yq(:,1);philq=Yq(:,2);Zlq=Yq(:,3);

    % Find the corresponding indexes of the extreme max values 
    extrMaxIndex=find(diff(sign(diff(Zlq)))==-2)+1;
    m=length(extrMaxIndex)-1; dtheta=m*2*pi;
    philq1=philq(extrMaxIndex(1));philq2=philq(extrMaxIndex(end));
    dphi=philq2-philq1;
    qq(iq)=dphi/dtheta;
    psiq(iq)=fPsi(Rq,Zq);
end
%     qq
%     psiq
subplot(224); plot(psiq,qq,'LineWidth',2); grid on;
ylim([0,1.2*max(qq)]); % xlim([0,max(psiq)]);
title(['q(\psi), \Psi_0=',num2str(Psi0)]);

print('-dpng',['SolovevEquilibrium,kappa=',num2str(kappa),...
    ',Psi0=',num2str(Psi0),',Ra0=',num2str(Ra0/R0),...
    ',B0=',num2str(B0),'.png']);
end

function dy=fieldline(lp,y) % cal field line trajectory
    global fBR fBZ fBp fBt;
    dy=zeros(3,1); % R phi Z
    Rl=y(1);phil=y(2);Zl=y(3);
    dy=[fBR(Rl,Zl)/fBp(Rl,Zl);
        fBt(Rl,Zl)/(fBp(Rl,Zl)*Rl);
        fBZ(Rl,Zl)/fBp(Rl,Zl)];    
end