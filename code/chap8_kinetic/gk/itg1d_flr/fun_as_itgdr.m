% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2016-05-09 15:10
% Adaptive Simpson quadrature formula to calculate the 2D integral
% 16-09-30 08:30 update
function Inm=fun_as_itgdr(w,ky,kapn,kapt,ws,wd)

% n=1; m=0; w=0.0002+0.001i; k=5;

% fnm=@(x,y) 2/sqrt(pi)*besselj(0,k*abs(y)).^2.*x.^m.*abs(...
%     y).^n./(x.^2+y.^2/2+w).*exp(-(x.^2+y.^2));
fnm=@(x,y) 2/sqrt(pi)*besselj(0,sqrt(2)*ky*abs(y)).^2.*(w-ws*(kapn+...
    kapt*(x.^2+y.^2-3/2)))./(w-2*wd*(x.^2+y.^2/2)).*abs(y).*exp(-(x.^2+y.^2));

% xmax=3e1; xmin=-0.0;  ymin=0; ymax=3e1; tol=1e-8;
% Inm=2*dblquad(fnm,xmin,xmax,ymin,ymax,tol);

xmax=1.2e1; xmin=-xmax;  ymin=0; ymax=2e1; tol=1e-6;
Inm=dblquad(fnm,xmin,xmax,ymin,ymax,tol);
