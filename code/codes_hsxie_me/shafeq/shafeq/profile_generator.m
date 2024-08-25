% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-03-10 08:14
% Tokamak q, n, T profile generator, for Shafranov shift equilibrium
% Normalization: White2001 book, B0=1, R0=1

close all; clear; clc;
opt=3;

% % eq=8, xhs ballooning
% r0=66.6; % major radius, unit=cm
% b0=16772; % on-axis magnetic field, unit=gauss
% etemp0=1.5625e3; % on-axis electron temperature, unit=ev
% % eden0=1.5346e14/3; % on-axis electron number density, unit=1/cm^3
% eden0=0.5e14;

% % eq=8, xhs ballooning
r0=165; % major radius, unit=cm
b0=13500; % on-axis magnetic field, unit=gauss
etemp0=0.4e3; % on-axis electron temperature, unit=ev
% eden0=1.5346e14/3; % on-axis electron number density, unit=1/cm^3
eden0=5e13;

% % eq=1, Cyclone
% r0=83.5;
% b0=20125.4;
% etemp0=2223.0;
% eden0=0.3142e14;

% % eq=3, Wenjun RSAE case
% r0=141.94; % major radius, unit=cm
% b0=19100.0; % on-axis magnetic field, unit=gauss
% etemp0=2500.0; % on-axis electron temperature, unit=ev
% eden0=1.448e14/2;

% etemp0=etemp0/2; eden0=eden0/2;

itemp0=etemp0;
iden0=eden0;
ftemp0=2.0*etemp0;
% fden0=1.0e-5*eden0;
fden0=1.0e-3*eden0;

% keep this the same with GTC analytic.F90 & eqdata.F90
fload=0; qelectron=-1; qion=1; qfast=1;
if(fload==0)
    nipp=abs(qelectron)/qion;
    fden0=1.0e-6*eden0; % eqdata.F90, subroutine analyticeq
else
    nipp=(abs(qelectron)-fden0*qfast)/qion;
end

%%
% npp=50; 
npp=256;
if(opt==1) 
    % 1, original GTC analytical eq, q=q1+q2*psi+q3*psi^2, n,T=tanh()

% % eq=1, Cyclone    
%     psiw=0.0375;
%     q1=0.82; q2=1.1; q3=1.0;
%     ne1=0.205; ne2=0.3; ne3=0.4;
%     te1=0.415; te2=0.18; te3=0.4;

% % eq=3, Wenjun RSAE case
%     psiw=0.02895;
%     q1=1.87; q2=-1.2; q3=2.0;
%     ne1=0.0; ne2=0.0; ne3=1;
%     te1=0.0; te2=0.0; te3=1.0;

% % eq=8, xhs ballooning
%     psiw=0.012;
%     q1=1.5; q2=0.0; q3=1.0;
%     psiw=0.0075;
%     q1=1.2;q2=1.5;q3=0.9;
    psiw=0.012;
    q1=1.0; q2=0.0; q3=2.0;
    ne1=0.4944; ne2=0.7; ne3=0.04;
    te1=0.0; te2=0.0; te3=1.0;
    
%     psiw=0.05; % psi on wall
    psi=0:(1/npp)*psiw:psiw;
    psi_n=psi./psiw; % normalization psi

%     q1=1.0; q2=0.0; q3=2.0;
    q=q1+q2.*psi_n+q3.*psi_n.^2;

    r=sqrt(2.*(q1.*psi_n+(q2.*psi_n.^2)./2.0+(q3.*psi_n.^3)./3.0)*psiw);
    a=max(r); % minor radius, a/R0
    r_n=r./a; % r/a

%     ne1=0.4944; ne2=0.7; ne3=0.05;
    ne_n=1.0+ne1.*(tanh((ne2-(psi_n))./ne3)-1.0);

%     te1=0.0; te2=0.0; te3=1.0;
    te_n=1.0+te1.*(tanh((te2-(psi_n))./te3)-1.0);
    
    ti_n=te_n;
    nf_n=ne_n;
    tf_n=ti_n;
    ni_n=ne_n;

elseif(opt==2 || opt==3)
    r_n=0:(1/npp):1;
    dr=r_n(2)-r_n(1);
    
    if(opt==2)
        % 2, 11-poits fitting q, n, T
% qdat=[0  ,  0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ;
%       1.00, 1.05, 1.13, 1.35 , 1.76, 2.15 , 2.67, 3.25 , 3.82, 4.45, 5.02];
%   
%         q=interp1(qdat(1,:),qdat(2,:),r_n,'spline');

%         q0=1.4; q1=1.9; q2=1.2; q3=3-(q0+q1+q2);
%         q=q0+q1*r_n+q2*r_n.^2+q3*r_n.^2;
        
        q=1.2+1.3*(tanh((r_n-0.5)/1.2)+1);
        
  
    else
        % 3, 11-poits fitting s, ...
% sdat=[0  ,  0.2 , 0.4 , 0.6 , 0.75 , 0.81, 0.88 , 0.91 , 0.94 , 0.97, 1.0 ;
%       0.00, 0.1,  0.3,  0.4 , 0.5,  0.61 , 0.81, 0.91 , 0.95, 0.98, 1.0];
%         s=interp1(sdat(1,:),sdat(2,:),r_n,'cubic');
        s=1-1*(r_n-1).^2;

        q0=0.7; % set a guess q(r=0)
        tmp=[0,cumsum(s(2:end)./r_n(2:end)).*dr];
        q=q0.*exp(tmp);
        
%         qa=3.5; % fix qa
        qa=3.0; % fix qa
        q=q./q(end)*qa;
        
    end
    
    psi_tmp=[0,cumsum(r_n(2:end)./q(2:end)).*dr];  % psip=integral(r/q)
    a=40/165; % set a=r0/R0
    psi=a^2*psi_tmp;
    psiw=max(psi);
    psi_n=psi./psiw; % normalization psi
    r=r_n.*a;
    
    % 11-poits fitting n, T
% nedat=[0  ,  0.2 , 0.4 , 0.6 , 0.72 , 0.81, 0.85 , 0.875 , 0.905 , 0.93, 0.97, 1.0 ;
%       1.00, 0.99, 0.97, 0.93 , 0.89 , 0.85, 0.75 , 0.51, 0.21 , 0.09, 0.045, 0.04];
% nedat=[0  ,  0.2 , 0.4 , 0.56 , 0.68 , 0.81, 0.85 , 0.875 , 0.905 , 0.93, 0.97, 1.0 ;
%       1.00, 0.9, 0.7, 0.5 , 0.4 , 0.3, 0.2 , 0.15, 0.105 , 0.08, 0.055, 0.05]; % 2014-07-27
nedat=[0  ,  0.2 , 0.4 , 0.56 , 0.68 , 0.84, 0.87 , 0.90 , 0.925 , 0.945, 0.965,  1.0 ;
      1.00,  0.9,  0.79,  0.725 , 0.65, 0.55, 0.35 , 0.21,  0.135 ,  0.10, 0.065, 0.05];
% tedat=[0  ,  0.2 , 0.4 , 0.6 , 0.75 , 0.88, 0.91 , 0.93 , 0.955 , 0.98, 1.0 ;
%       1.00, 0.99, 0.97, 0.93 , 0.89 , 0.85, 0.79 , 0.55, 0.29 , 0.15, 0.1];
    tedat=nedat;
%   tedat(2,:)=0.*tedat(2,:)+1.0;
   
%     ne_n=interp1(nedat(1,:),nedat(2,:),r_n,'spline');
%     te_n=interp1(tedat(1,:),tedat(2,:),r_n,'spline');
    ne_n=interp1(nedat(1,:),nedat(2,:),r_n,'cubic');
    te_n=interp1(tedat(1,:),tedat(2,:),r_n,'cubic');
    
    ni_n=ne_n;
    ti_n=te_n;
    nf_n=ne_n;
    tf_n=te_n;
    
else
    % Reading data from TRANSP
    % copy and modify from previous file "nc_readncwp2.m"
    % need rewrite
    
    ncid=netcdf.open('iterdb.nc.01550','NC_NOWRITE');
    %load profiletmp;
    load -ascii proshape1.dat;

    %mpsi=length(psir);
    dimid=netcdf.inqDimID(ncid,'dim_rho');
    [dimname,mpsi]=netcdf.inqDim(ncid,dimid);

    lsp=100; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    len2=15;

    rad=[1:mpsi]';
    rad1=linspace(1,mpsi,lsp)';

    dataid=netcdf.inqVarID(ncid,'psir_grid');%%%%%%%%%%
    psir=netcdf.getVar(ncid,dataid);%%%%%%%%%%%%%%%

    % psipol1=interp1(rad,psir,rad1,'spline');
    psiptmp=linspace(psir(1),psir(end),lsp+1);
    psipol1=psiptmp(2:end);

    dataid=netcdf.inqVarID(ncid,'rho_grid');%%%%%%%%%
    rho_grid=netcdf.getVar(ncid,dataid);%%%%%%%%%%%

    dataid=netcdf.inqVarID(ncid,'Te');%%%%%%%%%
    Te=netcdf.getVar(ncid,dataid);%%%%%%%%%%%

    dataid=netcdf.inqVarID(ncid,'ene');%%%%%%%%%
    ene=netcdf.getVar(ncid,dataid);%%%%%%%%%%%

    dataid=netcdf.inqVarID(ncid,'Ti');%%%%%%%%%
    Ti=netcdf.getVar(ncid,dataid);%%%%%%%%%%%

    dataid=netcdf.inqVarID(ncid,'enion');%%%%%%%%%
    enitemp=netcdf.getVar(ncid,dataid);%%%%%%%%%%%
    eni=enitemp(:,1);%?????????????????????????????????????????
    enimp=enitemp(:,2);%?????????????????????????????????????????
    %??????????????????????????????????????????????????????????????????
    %dataid=netcdf.inqVarID(ncid,'tglf_imp_e_flux');%%%%%%%%%
    %enimp=netcdf.getVar(ncid,dataid);%%%%%%%%%%%
    %??????????????????????????????????????????????????????????????????

    profile=zeros(mpsi,len2);
    profile1=zeros(lsp,len2);
    profile1(:,1)=psipol1;

    profile(:,1)=psir;
    profile(:,2)=rho_grid;
    profile(:,3)=rho_grid;
    profile(:,4)=1.0;
    profile(:,5)=1.0+rho_grid;
    profile(:,6)=Te;
    profile(:,7)=ene;
    profile(:,8)=Ti;
    profile(:,9)=1.0;
    profile(:,10)=0.0;
    profile(:,11)=0.0;
    profile(:,12)=eni;
    profile(:,13)=enimp;%??????????????????????????????
    profile(:,14)=1.0;
    profile(:,15)=1.0;

    %estimate radial domain size
    rc=rho_grid(end)*0.8; %in meter
    psic=interp1(rho_grid,psir,rc);
    rho1=interp1(psir,rho_grid,0.27);
    rho0=2*rc-rho1;
    psi00=interp1(rho_grid,psir,rho0);
    psi11=interp1(rho_grid,psir,rho1);

    len1=length(psir);
    indx1=sum(psipol1<psi00);
    indxc=sum(psipol1<psic);
    indx2=len1-sum(psipol1>psi11);

    psi0=interp1(rho_grid,psir,rho0-0.02);
    psi0/psir(end);
    psi1=interp1(rho_grid,psir,rho1+0.02);
    psi1/psir(end);

    tec=interp1(rho_grid,Te,rc);
    tic=interp1(rho_grid,Ti,rc);
    nec=interp1(rho_grid,ene,rc);

    for j=2:len2,
        profile1(:,j)=interp1(profile(:,1),profile(:,j),psipol1,'spline');
    end

    %smooth function at two ends
    dpsi=psipol1(2)-psipol1(1);

    for jt=6:8
    profiletmp=profile1(:,jt);
    profileint=zeros(size(profiletmp));
    for j=2:length(profiletmp)
        profileint(j)=profileint(j-1)+(proshape1(j-1,2)*profiletmp(j-1)+proshape1(j-1,2)*profiletmp(j-1))*0.5*dpsi;
    end
    profiletmp1=profiletmp.*proshape1(:,2)-profileint;
    profiletmp1=profiletmp1(indxc)-profile1(indxc,jt);
    profile1(:,jt)=profiletmp1;
    end

    %%correction to ni to make ni=ne
    %%normalize psipol1
    profile1(:,1)=profile1(:,1)/profile1(end,1);
    profile1(:,12)=profile1(:,7);
    
    
end

if(opt==1 || opt==2 || opt==3)
    lsp=length(r);

    Zeff=zeros(1,lsp);
    Er=zeros(1,lsp);
    nimp=zeros(1,lsp);
    omega_tor=zeros(1,lsp);

    x=sqrt(psi);

    % s=(r/q)*(dq/dr)
    s=r_n(1:end-1).*(log(q(2:end))-log(q(1:end-1)))./diff(r_n);
%     s=[s(1),s];
    s=[s,2*s(end)-s(end-1)];

    R0=r0; % unit=cm
%     r=r_n*r0; % unit=cm
    r=r_n*a*R0; % unit=cm
    R=R0+r; % unit=cm
    ne=ne_n*eden0; % unit=1/cm^3
    te=te_n*etemp0; % unit=ev
    ni=ni_n*iden0;
    ti=ti_n*itemp0;
    nf=nf_n*fden0;
    tf=tf_n*ftemp0;
    
    rho_M=ni/ni(1); % mass density

    betae=4.03e-11.*te.*ne/(b0*b0); % average <bb> is equal to b0=1
    betai=4.03e-11.*ti.*ni/(b0*b0);
    betaf=4.03e-11.*tf.*nf/(b0*b0);
%     beta=betae+betai+betaf;
    beta=betae;
    p=beta;

    % alpha=-q^2*R*(dbeta/dr), R=1
    alpha=-R0*q(1:end-1).^2.*diff(beta)./diff(r);
    alpha=[alpha(1),alpha];

end

%% plot
hf=figure('unit','normalized','Position',[0.1 0.4 0.7 0.53],...
            'Name','Equilibrium profile',... % 'menubar','none',...
            'NumberTitle','off');
set(gcf,'DefaultAxesFontSize',14);

subplot(341); plot(psi_n,q,'--g',r_n,q,'r','LineWidth',2);
xlim([0,1]); ylim([0,5]); grid on;
title('q');
legend('(psi_n)','(r_n)',2); legend('boxoff');

subplot(342);plot(psi_n,r,'--g','LineWidth',2);
grid on;xlim([0,1]);title('r(psi_n)');

subplot(343); plot(psi_n,ne_n,'--g',r_n,ne_n,'r','LineWidth',2);
xlim([0,1]); ylim([0,1.2]); grid minor;
title('ne');
subplot(344);plot(psi_n,te_n,'--g',r_n,te_n,'r','LineWidth',2);
xlim([0,1]);grid on;
title('Te');

subplot(345); plot(psi_n,ni_n,'--g',r_n,ni_n,'r','LineWidth',2);
xlim([0,1]); ylim([0,1.2]); grid minor;
title('ni');
subplot(346);plot(psi_n,ti_n,'--g',r_n,ti_n,'r','LineWidth',2);
xlim([0,1]);grid on;
title('Ti');

subplot(347); plot(psi_n,nf_n,'--g',r_n,nf_n,'r','LineWidth',2);
xlim([0,1]); ylim([0,1.2]); grid minor;
title('nf');
subplot(348);plot(psi_n,tf_n,'--g',r_n,tf_n,'r','LineWidth',2);
xlim([0,1]);grid on;
title('Tf');

subplot(3,4,9);
plot(psi_n,s,'--g',r_n,s,'r','LineWidth',2);
xlim([0,1]); title('s=(r/q)(dq/dr)'); grid minor;

subplot(3,4,10);
plot(psi_n,alpha,'--g',r_n,alpha,'r','LineWidth',2);
xlim([0,1]);title('\alpha=-q^2R(d\beta/dr)');grid minor;

subplot(3,4,11);
dlnne=-(log(ne(2:end))-log(ne(1:end-1)))./(r_n(2:end)-r_n(1:end-1));
plot(psi_n(1:end-1),-(log(ne(2:end))-log(ne(1:end-1)))./(psi_n(2:end)-psi_n(1:end-1)),...
    '--g',r_n(1:end-1),dlnne,'r','LineWidth',2);
xlim([0,1]);title('-dln(ne)');grid on;

print('-dpng','profile_dat.png');


%% output

fid_tmp = fopen('profile_tmp.dat', 'w'); % for shafeq_spdata.m
profile_tmp=[R; r; q; p];
fprintf(fid_tmp,'%16.8e %16.8e %16.8e %16.8e\n',profile_tmp);
fclose(fid_tmp);

awcon=1;
if(awcon==1)
    fid_sa = fopen('profile_sa.dat', 'w'); % for awcon.m
    profile_sa=[r_n; q; p; rho_M];
    fprintf(fid_sa,'%16.8e\n',a);
    fprintf(fid_sa,'%16.8e %16.8e %16.8e %16.8e\n',profile_sa);
    fclose(fid_sa);
end

% type profile_tmp.dat
zeff=1+0.*x; Er=0+0.*x; omega_tor=0+0.*x; nimp=0+0.*x;
% nf=0+0.*x; tf=0+0.*x;
fid = fopen('profile.dat', 'w'); % for GTC
profile=[psi_n; x; r; R0+0.*x; R; te; ne; ti; zeff; omega_tor; Er; ni; nimp; nf; tf];
fprintf(fid,'%s\n','    Pol-Flux             x               r               R                R+r              Te                ne               Ti             Zeff           omega-tor            Er               ni             nimp              nf               Tf');
fprintf(fid,'%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n',profile);

strend=[10,10,10,'Poloidal Flux in Webers',10,...
'x is square root of normalized toroidal flux',10,...
'r is midplane radius in cm (midplane halfwidth of flux surface width at midplane)',10,...
'R is flux surface center in cm',10,...
'R+r is local major radius in cm on the outer midplane',10,...
'Te is electron temperature in eV',10,...
'ne is electron density in m^-3',10,...
'Ti is ion temperature in eV; last two points are artificial',10,...
'Zeff is from Carbon density profile measurement',10,...
'omega-tor is measured angular velocity from carbon rotation in radians/sec;to get votr, multiply omega by local R from equilibrium (or from R+r)',10,...
'Er is radial electric field in volts/cm',10,...
'ni is ion density in m^-3',10,...
'nimp is impurity density in 10^19 m^-3',10,...
'nf is fast ion density in m^-3',10,...
'Tf is electron temperature in eV'];

fprintf(fid,'%s\n',strend);


fclose(fid);

