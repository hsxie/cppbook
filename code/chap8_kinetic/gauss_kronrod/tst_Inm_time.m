% 16-05-23 20:34
close all; clear; clc;
format long;

n=1; m=2; b=0.8; kpar=0.1;
w=0.002+0.004i; 
% w=-0.002-0.004i; 

nt1=100;
time1=cputime;
for jt=1:nt1
    fInm=fun_gz_gk_Inm(w,kpar,b,n,m);
end
tgl=(cputime-time1)/nt1
fInm

fnm=@(x,y) 2/sqrt(pi)*besselj(0,sqrt(2*b)*abs(y)).^2.*x.^m.*abs(y).^n./(x.^2+...
    y.^2/2+w-kpar*x).*exp(-(x.^2+y.^2));

xmax=1e1; xmin=-xmax;  ymin=0; ymax=2e1; tol=1e-8;
nt2=2;
time3=cputime;
for jt=1:nt2
    Inm=dblquad(fnm,xmin,xmax,ymin,ymax,tol);
end
tas=(cputime-time3)/nt2
Inm
%%
% Inm/fInm