% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2016-10-10 11:57
% Adaptive Simpson quadrature formula to calculate the 2D integral
% for entropy mode.
% 16-10-16 16:14 with kpara
% 16-10-18 23:24 fixed a bug for kpara \neq 0
function fdr=fun_entropy_dr(w,k,wd,kapn,kapt,kz,tol,xmax,ymax)

if (nargin<9)
    ymax=1e1;
end
if (nargin<8)
    xmax=1e1;
end
if (nargin<7) % control the accuracy
    tol=1e-10;
end
if (nargin<6)
    kz=0.0;
end

nw=length(w); fdr=0.*w;

for jw=1:nw
	f=@(x,y) besselj(0,k*abs(y)).^2.*(w(jw)-wd*(kapn+...
		kapt*((x.^2+y.^2)/2-3/2))).*abs(y)./(w(jw)-kz*x-wd*(x.^2+...
		y.^2/2)).*exp(-(x.^2+y.^2)/2);

	xmin=-xmax;  ymin=0; 
	fdr(jw)=1.0-1.0/sqrt(2*pi)*dblquad(f,xmin,xmax,ymin,ymax,tol);
end

