clear;clc;
nii=5; njj=5; ni=2^nii+1; nj=2^njj;
dx=1.0; dz=2.0;
gamma=1.66667;
eta=0.01;
nu=0.05;
beta=0.5;
rho0=1.0;
b0=1.0;
bl=3.0;
va=sqrt(b0*b0/rho0);
p0=0.5*beta*b0*b0;
t0=0.5*beta*va*va;

x=zeros(ni,nj,5);
for i=1:ni
    s=(i-ni)*dx/bl;
    s1=b0*bl*log(cosh(s));
    b=b0*tanh(s);
    p=p0+0.5*(b0^2-b^2);
    rho=p/t0;
    for j=1:nj
        x(i,j,5)=s1;
        x(i,j,1)=rho;
        x(i,j,2)=p;
    end
end
rhoo=x(:,:,1);pp=x(:,:,2);uxx=x(:,:,3);uzz=x(:,:,4);ayy=x(:,:,5);
pcolor(ayy); shading('interp');