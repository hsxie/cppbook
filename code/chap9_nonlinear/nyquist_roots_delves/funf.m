% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-10-03 23:40
% funf.m, give f(z)=D'/D which determined by D(z)
function out=funf(z)

global mi Ti vti me tau Te vte alphai alphae theta k kappan kappat;

fnum=3;
if(fnum==0)
    % D=z^2*(z-1)
    out=(3*z.^2-2*z)./(z.^2.*(z-1));
elseif(fnum==1)
    % D=z^3-1
    out=(3.*z.^2)./(z.^3-1);
elseif(fnum==2)
    % D=cos(z)-z;
    out=(-sin(z)-1.0)./(cos(z)-z);
%     % D=sin(z)-z;
%     out=(cos(z)-1.0)./(sin(z)-z);
elseif(fnum==3) % Landau damping
    k=1.0/2;
    zeta=@(z)faddeeva(z,64)*1i*sqrt(pi);
    % f_fun=@(z)k*k+(1+(z/k/sqrt(2)).*zeta(z/k/sqrt(2)));
    out=((zeta(z/k/sqrt(2))-2*(z/k/sqrt(2)).*(1+...
        (z/k/sqrt(2)).*zeta(z/k/sqrt(2))))/k/sqrt(2))./(k*k+...
        (1+(z/k/sqrt(2)).*zeta(z/k/sqrt(2))));
elseif(fnum==4)
        
%     mi=1.0; Ti=1.0; vti=sqrt(Ti/mi);
%     me=1/1837; tau=1.0; Te=tau*Ti; vte=sqrt(Te/me);
%     alphai=1; alphae=-1/me; theta=0.01; 
%     k=1.0; kappan=0.0; kappat=2.0;

%     zeta=@(x)faddeeva(x)*1i*sqrt(pi);
    zeta=@(x)faddeeva(x,64)*1i*sqrt(pi);
    withdf=1;
    if(withdf==1)
       out=(1/k)*((-2*(kappan-kappat)*(1.0+z/(sqrt(2)*k*theta*vti).*zeta(z/(sqrt(2)*k*theta*vti)))+...
        (kappat*z/(sqrt(2)*k*theta*vti)+sqrt(2)*alphai*theta/vti).*(zeta(z/(sqrt(2)*k*theta*vti))-...
        2*z/(sqrt(2)*k*theta*vti).*(1+(z/(sqrt(2)*k*theta*vti)).*zeta(z/(sqrt(2)*k*theta*vti)))))/(sqrt(2)*theta*vti)^2-...
        (-2*(kappan-kappat)*(1.0+z/(sqrt(2)*k*theta*vte).*zeta(z/(sqrt(2)*k*theta*vte)))+...
        (kappat*z/(sqrt(2)*k*theta*vte)+sqrt(2)*alphae*theta/vte).*(zeta(z/(sqrt(2)*k*theta*vte))-...
        2*z/(sqrt(2)*k*theta*vte).*(1+(z/(sqrt(2)*k*theta*vte)).*zeta(z/(sqrt(2)*k*theta*vte)))))/(sqrt(2)*theta*vte)^2)./(...
        k*k+((kappan-kappat/2)*zeta(z/(sqrt(2)*k*theta*vti))+...
        (kappat*z/(sqrt(2)*k*theta*vti)+sqrt(2)*alphai*theta/vti).*(1+...
        z/(sqrt(2)*k*theta*vti).*zeta(z/(sqrt(2)*k*theta*vti))))/(sqrt(2)*theta*vti)-...
        ((kappan-kappat/2)*zeta(z/(sqrt(2)*k*theta*vte))+...
        (kappat*z/(sqrt(2)*k*theta*vte)+sqrt(2)*alphae*theta/vte).*(1+...
        z/(sqrt(2)*k*theta*vte).*zeta(z/(sqrt(2)*k*theta*vte))))/(sqrt(2)*theta*vte));
    else
        f=@(z)k*k+((kappan-kappat/2)*zeta(z/(sqrt(2)*k*theta*vti))+...
        (kappat*z/(sqrt(2)*k*theta*vti)+sqrt(2)*alphai*theta/vti).*(1+...
        z/(sqrt(2)*k*theta*vti).*zeta(z/(sqrt(2)*k*theta*vti))))/(sqrt(2)*theta*vti)-...
        ((kappan-kappat/2)*zeta(z/(sqrt(2)*k*theta*vte))+...
        (kappat*z/(sqrt(2)*k*theta*vte)+sqrt(2)*alphae*theta/vte).*(1+...
        z/(sqrt(2)*k*theta*vte).*zeta(z/(sqrt(2)*k*theta*vte))))/(sqrt(2)*theta*vte);
        eps=1e-6;
        out=(f(z+eps)-f(z))./f(z)/eps;
    end
end
