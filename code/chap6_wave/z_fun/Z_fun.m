% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2017-04-07 07:56
function fz=Z_fun(zeta)
[nx,ny]=size(zeta); fz=0.*zeta;
xmax=10; xmin=-xmax+0.1; tol=1e-5;
for jx=1:nx
    for jy=1:ny
        f=@(x) exp(-x.^2)./(x-zeta(jx,jy));
        fz(jx,jy)=1/sqrt(pi)*quad(f,xmin,xmax,tol);
        if(imag(zeta(jx,jy))<0)
            sigma=0;
%             sigma=2;
            fz(jx,jy)=fz(jx,jy)+1i*sigma*sqrt(pi)*exp(-zeta(jx,jy)^2);
        end
    end
end