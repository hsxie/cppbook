% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-03-10 10:24
% Shafranov shift equilibrium for Tokamak, 'profile.dat' -> 'spdata.dat'
% Normalization: White2001 book, B0=1, R0=1
% Ref/For: ORBIT, GTC, ... Meiss1990
% 2014-04-03 12:49 update & checked with GTC

% prepare "profile_tmp.dat" using "profile_generator.m" 
% "profile_tmp.dat" is for "shafeq_spdata.m" while "profile.dat" is for GTC
% You can choose your own other profiles.

% good approximation when 'eps=a/R < 0.25' & 'beta ~< eps^2'

% run profile_generator;

close all; clear; clc;

% ishift==1, with Shafranov shift; 0, not shift (GTC old analyticeq)
% =2, when p'=0, make shift=0 to model s-alpha, 2014-12-07 21:17
% =3, 1/3 of ishift=2, 2014-12-08 12:05
ishift=0; 

% ithetaf=0, old theta_s -> theta_flux, ishift=0 is ok, but 
% bad for ishift\=0; else, interp, more accurate
ithetaf=1;
fm=1; % fieldmodel=0 or 1 in GTC, 2014-11-06 09:22

% profile_tmp=importdata('profile_tmp.dat',' ',2); % read "profile_tmp.dat"
profile_tmp=importdata('profile_tmp.dat',' ');
R0=profile_tmp(1,1); % major radius, unit=cm
r_tmp=profile_tmp(:,2); % minor radius, unit=cm
q_tmp=profile_tmp(:,3); % q
p_tmp=profile_tmp(:,4); % beta

a=max(r_tmp)/R0; % minor radius, a/R0
r=0.0:0.001*a:a;
dr=r(2)-r(1);
nr=length(r);

r_tmp=r_tmp/R0;
q=interp1(r_tmp,q_tmp,r,'cubic');
p=interp1(r_tmp,p_tmp,r,'spline');

psi=[0,cumsum(r(2:end)./q(2:end)).*dr];  % psip=integral(r/q)

psiw=max(psi);
psi_n=psi./psiw; % normalization psi
r_n=r./a; % r/a

if(ishift~=0)
% cal the Shafranov shift, Delta' (which should < 1, since Delta'*dr<dr to 
% aviod the intersecting of two fluxs) and Delta
    sftmp=(ishift==1); % 2014-12-07 21:17
    Deltap=0.*r; Delta=0.*r;
    inttmp=0.*r; betap=0.*r; li=0.*r;
    % an analytical expression is not easy, so we use numerical integral
    for ir=3:nr-1 % singularity at r=0
        dr=r(ir)-r(ir-1);
    %     rtmp=(r(ir)+r(ir-1))/2; qtmp=(q(ir)+q(ir-1))/2;
        rtmp=r(ir); qtmp=q(ir);
        pptmp=(p(ir+1)-p(ir-1))/dr/2;
%         inttmp(ir)=inttmp(ir-1)+(rtmp^2/qtmp^2-2*rtmp*pptmp)*rtmp*dr;
        % R. B. White book, p=beta/2
        inttmp(ir)=inttmp(ir-1)+(sftmp*rtmp^2/qtmp^2-1*rtmp*pptmp)*rtmp*dr;
    %     inttmp(ir)=inttmp(ir-1)+(-2*rtmp*pptmp)*rtmp*dr;
    %     inttmp(ir)=inttmp(ir-1)+(rtmp^2/qtmp^2)*rtmp*dr;
%         Deltap(ir)=inttmp(ir)*qtmp^2/rtmp^3;
        Deltap(ir)=inttmp(ir)*qtmp^2/rtmp^3;
        Delta(ir)=Delta(ir-1)+Deltap(ir)*dr;
    end
    Deltap(nr)=Deltap(nr-1);
    Delta(nr)=Delta(nr-1);
    
    if(ishift>=3) % 2014-12-08 12:05
        Deltap=Deltap/ishift;
        Delta=Delta/ishift;
    end

    inttmp=0.*r;
    for ir=nr-1:-1:2
        dr=r(ir+1)-r(ir);
        qtmp=q(ir);
        r2qptmp=(r(ir+1)^2/q(ir+1)-r(ir-1)^2/q(ir-1))/dr/2;
        inttmp(ir)=inttmp(ir+1)+(r2qptmp/qtmp)*dr;
        g2(ir)=-p(ir)+inttmp(ir);
    end
    g2(1)=g2(2); g2(nr)=g2(nr-1);
else
    Deltap=0.*r; 
    Delta=0.*r;
    g2=0.*r;
end
g=1+g2;    


%% plot
hf=figure('unit','normalized','Position',[0.01 0.5 0.5 0.43],...
            'Name','Equilibrium profile',... % 'menubar','none',...
            'NumberTitle','off');
set(gcf,'DefaultAxesFontSize',14);

q55=q; r55=r;

subplot(231); plot(psi_n,q,'--g',r_n,q,'r','LineWidth',2);
xlim([0,1]); ylim([0,5]); grid on;
title('q');
legend('(psi_n)','(r_n)',2); legend('boxoff');

subplot(232);plot(psi_n,r,'--g','LineWidth',2);
% grid on;
xlim([0,1]);title('r(psi_n)');

subplot(233); plot(psi_n,p,'--g',r_n,p,'r','LineWidth',2);
xlim([0,1]); ylim([0,max(p)]); 
% grid minor;
title('p');

subplot(234); plot(psi_n,g,'--g',r_n,g,'r','LineWidth',2);
xlim([0,1]); ylim([0.8,1.2]); 
% grid minor;
title('g=1+g2');

subplot(235); plot(psi_n,Delta,'r',psi_n,Deltap,'--g',[0,1],[1,1],'--','LineWidth',2);
% grid on;
xlim([0,1]); axis tight;
title('\Delta(psi_n)');
legend('\Delta','\Delta''',2); legend('boxoff');

yy1=Delta; yy2=Deltap;
ydat2=[r_n;yy1;yy2]; save ydat2;

%% cal f1, ..., f9, ...
% dthe=pi/50;
dthe=pi/250;
the=dthe:dthe:2*pi;
nthe=length(the);
nr=length(r);
[rtmp,ttmp]=meshgrid(r,the);
ptmp=repmat(psi,nthe,1);
etatmp=(repmat(Deltap,nthe,1)+rtmp)./2;

if(ithetaf==0) % old method, ishift=0 is enough, order of epsilon^2
    % GTC eqdata.F90, subroutine analyticeq:
    %          B=1-rcos(theta_geo)+(rcos(theta_geo))**2 to order of epsilon^2

    % Xtmp=1+rtmp.*cos(ttmp)-repmat(Delta,nthe,1);
    % Ztmp=rtmp.*sin(ttmp);

    % below, ttmp is theta_flux, 'theta_geo -> theta_flux' has included
    Xtmp=1+rtmp.*cos(ttmp)-repmat(Delta,nthe,1)+fm*etatmp.*rtmp.*(cos(2*ttmp)-1);
    Ztmp=rtmp.*sin(ttmp)+fm*etatmp.*rtmp.*sin(2*ttmp);

    % Bptmp=rtmp./(repmat(q,nthe,1).*Xtmp.*(1-repmat(Deltap,nthe,1).*cos(ttmp)));
    % Bttmp=repmat(g,nthe,1)./Xtmp;
    % Btmp=sqrt(Bptmp.^2+Bttmp.^2);
    % Btmp=1-rtmp.*cos(ttmp)+(rtmp.*cos(ttmp)).^2;
    Btmp=1-rtmp.*cos(ttmp)-fm*etatmp.*rtmp.*(cos(2*ttmp)-1)+(rtmp.*cos(ttmp)).^2;

    Itmp=rtmp.^2./(repmat(q,nthe,1).*Xtmp.*(1-repmat(Deltap,nthe,1).*cos(ttmp)));
    % Itmp=rtmp.^2./(repmat(q,nthe,1)); % only first order

    Jactmp=repmat(q,nthe,1).*Xtmp.*(1-repmat(Deltap,nthe,1).*cos(ttmp));
else
    % interp method,  to avoid deformation of flux surface
    %  Boozer coordinate, solve, theta_f=theta_s-2*eta*sin(theta_s), 
    % or, theta_s=theta_f+2*eta*sin(theta_f)
    intp=1;
    if(intp==1)
        tstmp=[];
        for ir=1:nr
            thetaf=the-(r(ir)+Deltap(ir))*sin(the); % Boozer coordinate
            thetaftmp=the;
            % cal theta_s using interp, since the inverse function is not explicit
            theta=interp1(thetaf,the,thetaftmp);
            tstmp(:,ir)=theta;
        end
    else
        tstmp=ttmp+(repmat(Deltap,nthe,1)+rtmp).*sin(ttmp); % theta_s
    end
    Xtmp=1+rtmp.*cos(tstmp)-repmat(Delta,nthe,1);
    Ztmp=rtmp.*sin(tstmp);
    
    % below need check
    Btmp=repmat(g,nthe,1)./Xtmp;
%     Btmp=1-rtmp.*cos(ttmp)-etatmp.*rtmp.*(cos(2*ttmp)-1)+(rtmp.*cos(ttmp)).^2;
    Itmp=rtmp.^2./(repmat(q,nthe,1)); % only first order

    Jactmp=repmat(q,nthe,1).*Xtmp.*(1-repmat(Deltap,nthe,1).*cos(tstmp));
end
nutmp=-repmat(q,nthe,1).*(rtmp+repmat(Deltap,nthe,1)).*cos(ttmp); % need check
deltmp=0.*Xtmp;

krip=0; nrip=0; brip=0;  % ripple parameter, we set them zero here
dum1=0; % dum1 = wrip*rmaj/xc
dum2=0; % dum2 = xrip*rmaj/xc
d0=0;
hatmp=0.*Xtmp;
hbtmp=0.*Xtmp;

% set output parameters
rmaj=R0;   % magnetic axis in cm
lsp=81; % b,x,z,giac,q,ripple- poloidal spline grid points
% lst=81; % b,x,z,giac,ripple, must be odd-theta spline grid points
lst=161;
lemax=4; % g, I
lrmax=8; % pol(r), r(pol), pressure
ped=psiw;
pw=psiw;

% rt=linspace(0,a,lsp);
pt=linspace(0,psiw,lsp);
% pst=linspace(0,pw,lsp);
dthe=2*pi/(lst-0);
the=dthe:dthe:2*pi;
% the=linspace(0,2*pi,lst);
% [rr,tt]=meshgrid(rt,the);
[pp,tt]=meshgrid(pt,the);


X=interp2(ptmp,ttmp,Xtmp,pp,tt);
Z=interp2(ptmp,ttmp,Ztmp,pp,tt);
% X=interp2(rtmp,ttmp,Xtmp,rr,tt);
% Z=interp2(rtmp,ttmp,Ztmp,rr,tt);
% eta=(rr)./2;
% X=X+eta.*rr.*(cos(2*tt)-1);
% Z=Z+eta.*rr.*sin(2*tt);


I=interp2(ptmp,ttmp,Itmp,pp,tt);
B=interp2(ptmp,ttmp,Btmp,pp,tt);
Jac=interp2(ptmp,ttmp,Jactmp,pp,tt);
nu=interp2(ptmp,ttmp,nutmp,pp,tt);
del=interp2(ptmp,ttmp,deltmp,pp,tt);
ha=interp2(ptmp,ttmp,hatmp,pp,tt);
hb=interp2(ptmp,ttmp,hbtmp,pp,tt);

subplot(236);
surf(X,Z,I); xlabel('X'); ylabel('Z'); title('I');

print('-dpng','shafeq_spdat.png');

psit=r.^2/2;

q=interp1(psi,q,pt);
g=interp1(psi,g,pt,'spline');
p=interp1(psi,p,pt,'spline');
r=interp1(psi,r,pt,'spline'); %
psit=interp1(psi,psit,pt,'spline'); % toroidal psi, not (poloidal) psi 
[qd1,qd2,qd3]=spln1d(pt,q);
[gd1,gd2,gd3]=spln1d(pt,g);
[pd1,pd2,pd3]=spln1d(pt,p);
[rp1,rp2,rp3]=spln1d(pt,r);
[ps1,ps2,ps3]=spln1d(pt,psit);
[rd1,rd2,rd3]=spln1d(pt,sqrt(psi));
% GTC analytic eq, I-current=c-current in boozer coordinates, rd not used
% [rd1,rd2,rd3]=spln1d(pt,0.*psi); % why not used in GTC?

% [qd1,qd2,qd3]=construct_spline0(rt,q);
% [gd1,gd2,gd3]=construct_spline0(rt,g);
% [pd1,pd2,pd3]=construct_spline0(rt,p);
% [ps1,ps2,ps3]=construct_spline0(rt,psi);
% [rp1,rp2,rp3]=construct_spline0(rt,r);
% [rd1,rd2,rd3]=construct_spline0(rt,sqrt(psi));

ispln=0;
[b1,b2,b3,b4,b5,b6,b7,b8,b9]=spln2d(ispln,pt,the,B);
[x1,x2,x3,x4,x5,x6,x7,x8,x9]=spln2d(ispln,pt,the,X);
[z1,z2,z3,z4,z5,z6,z7,z8,z9]=spln2d(ispln,pt,the,Z);
[g1,g2,g3,g4,g5,g6,g7,g8,g9]=spln2d(ispln,pt,the,Jac);
% [g1,g2,g3,g4,g5,g6,g7,g8,g9]=spln2d(ispln,pt,the,0.*Jac); % =0, not used in GTC?
% [rd1,rd2,rd3,rd4,rd5,rd6,rd7,rd8,rd9]=spln2d(ispln,pt,the,I);
[nu1,nu2,nu3,nu4,nu5,nu6,nu7,nu8,nu9]=spln2d(ispln,pt,the,nu);
[dl1,dl2,dl3,dl4,dl5,dl6,dl7,dl8,dl9]=spln2d(ispln,pt,the,del);
[ha1,ha2,ha3,ha4,ha5,ha6,ha7,ha8,ha9]=spln2d(ispln,pt,the,ha);
[hb1,hb2,hb3,hb4,hb5,hb6,hb7,hb8,hb9]=spln2d(ispln,pt,the,hb);
[r1,r2,r3,r4,r5,r6,r7,r8,r9]=spln2d(ispln,pt,the,0.*X);


%% output "spdata.dat"
fid = fopen('spdata.dat', 'w'); % for GTC
fprintf(fid,'ShafranovShiftCircle \t hsxie_v01\n');
fprintf(fid,'%4d%4d%4d%4d\n',lsp,lst,lemax,lrmax);
fmat3='%18.9e%18.9e\n';
fmat4='%18.9e%18.9e%18.9e\n';
fprintf(fid,fmat3,pw,ped);
fmat='%18.9e%18.9e%18.9e%18.9e\n';
for j=1:lsp
    fprintf(fid,fmat,b1(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,b2(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,b3(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,b4(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,b5(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,b6(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,b7(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,b8(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,b9(:,j)); fprintf(fid,'\n');
    
    fprintf(fid,fmat,x1(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,x2(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,x3(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,x4(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,x5(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,x6(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,x7(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,x8(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,x9(:,j)); fprintf(fid,'\n');
    
    fprintf(fid,fmat,z1(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,z2(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,z3(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,z4(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,z5(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,z6(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,z7(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,z8(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,z9(:,j)); fprintf(fid,'\n');
    
    fprintf(fid,fmat,g1(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,g2(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,g3(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,g4(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,g5(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,g6(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,g7(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,g8(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,g9(:,j)); fprintf(fid,'\n');
%     
%     fprintf(fid,fmat,rd1(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,rd2(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,rd3(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,rd4(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,rd5(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,rd6(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,rd7(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,rd8(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,rd9(:,j)); fprintf(fid,'\n');
%     
%     fprintf(fid,fmat,nu1(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,nu2(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,nu3(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,nu4(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,nu5(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,nu6(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,nu7(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,nu8(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,nu9(:,j)); fprintf(fid,'\n');
%     
%     fprintf(fid,fmat,dl1(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,dl2(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,dl3(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,dl4(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,dl5(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,dl6(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,dl7(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,dl8(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,dl9(:,j)); fprintf(fid,'\n');    
    
    fprintf(fid,fmat4,qd1(j),qd2(j),qd3(j)); % q
    fprintf(fid,fmat4,gd1(j),gd2(j),gd3(j)); % g-current
    fprintf(fid,fmat4,rd1(j),rd2(j),rd3(j)); % i-current
    fprintf(fid,fmat4,pd1(j),pd2(j),pd3(j)); % pressure
    fprintf(fid,fmat4,rp1(j),rp2(j),rp3(j)); % minor radius
    fprintf(fid,fmat4,ps1(j),ps2(j),ps3(j)); % toroidal flux
end

fprintf(fid,'%4d%4d\n',krip,nrip);
fprintf(fid,'%18.9e%18.9e%18.9e\n',rmaj,d0,brip);
fprintf(fid,'%18.9e%18.9e\n',dum1,dum2);
for j=1:lsp
    fprintf(fid,fmat,r1(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,r2(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,r3(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,r4(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,r5(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,r6(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,r7(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,r8(:,j)); fprintf(fid,'\n');
    fprintf(fid,fmat,r9(:,j)); fprintf(fid,'\n');
    
%     fprintf(fid,fmat,ha1(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,ha2(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,ha3(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,ha4(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,ha5(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,ha6(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,ha7(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,ha8(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,ha9(:,j)); fprintf(fid,'\n');
%     
%     fprintf(fid,fmat,hb1(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,hb2(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,hb3(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,hb4(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,hb5(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,hb6(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,hb7(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,hb8(:,j)); fprintf(fid,'\n');
%     fprintf(fid,fmat,hb9(:,j)); fprintf(fid,'\n');
end

figure;
plot(X,Z,'.'); axis equal;
print('-dpng',['flux_surf_shift=',num2str(ishift),'.png']);

fclose(fid);

