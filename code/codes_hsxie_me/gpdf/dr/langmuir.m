% Hua-sheng XIE, 2011-10-08
% Dispersion relation solution for Langmiur wave Landau damping
clear;clc;
zeta=@(x)faddeeva(x)*1i*sqrt(pi);
f=@(x,k)1+k*k+x*zeta(x);w=[];
kmin=0.1;dk1=0.1;kmid=1;dk2=1;kmax=10.0;
k=[kmin:dk1:kmid,(kmid+dk2):dk2:kmax];
for kk=k
    options=optimset('Display','off');
    x=fsolve(f,1-0.1i,options,kk)*sqrt(2)*kk;
    w=[w,x];
end
wre=real(w);wie=imag(w);
wrt=1.0+1.5.*k.*k;
wit=-sqrt(pi/8).*exp(-1.0./(2.0.*k.^2)-1.5)./(k.^3);
loglog(k,wre,'-*r',k,-wie,'+r--',k,wrt,'b-',k,-wit,'b--');
legend('\omega_r','-\gamma','\omega_r','-\gamma','Location','SouthEast');grid on;
title(strcat('Langmuir Wave Dispersion Relation, approximate analytical'...
,10,'solution (dashed lines) and exact numerical computation (solidlines)'));
xlabel('k\lambda_D');ylabel('\omega/\omega_p');
xlim([kmin,kmax]);
ylim([0.0001,100]);