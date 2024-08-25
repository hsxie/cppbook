% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-10-03 22:39
% Using quad, instead of nx to calculate the integral contour
function sp=fun_sp(p,za,zb,tol)
% Numbers of zeros in rectangle domain, via Cauthy contour integral
% (xa,yb) ------------- (xb,yb)
%         |        *  |
%         |root(s)    |
%         |   *       |
% (xa,ya) ------------- (xb,ya)
% \int_c(f'/f)dz=2i*pi*N, N is the # of roots in complex domain.
% s_p=(1/2i*pi)*\int_c(z^p*f'/f)dz

if nargin<4, tol=[]; end
if isempty(tol), tol=1e-3; end

f=@(z)funf(z).*z.^p;
zc=real(zb)+1i*imag(za); zd=real(za)+1i*imag(zb);
sp=quad(f,za,zc,tol)+quad(f,zc,zb,tol)+quad(f,zb,zd,tol)+quad(f,zd,za,tol); 

sp=sp/(2*pi*1i);

