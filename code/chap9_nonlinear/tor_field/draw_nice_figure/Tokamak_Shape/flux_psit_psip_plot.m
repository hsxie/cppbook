% plot poloidal and toroidal flux psi_p and psi_t 
% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2013-01-12 08:02
close all; clear; clc;
figure; set(gcf,'DefaultAxesFontSize',15);

R0=1; a=0.2;

phimin=-0.0*pi; phimax=1.2*pi;
[Theta,Phi]=meshgrid(0:pi/20:2*pi,phimin:pi/20:phimax);
X31=(R0+a*cos(Theta+1.2*Phi)).*cos(Phi);
Y31=(R0+a*cos(Theta+1.2*Phi)).*sin(Phi);
Z31=a*sin(Theta+1.2*Phi);
F31=0.*X31+0.1;
surf(X31,Y31,Z31,F31,'EdgeColor','b','FaceColor','none');
view(22,28);axis equal; % plot the shape

% plot coordinates axis
arrow([0 0 -a],[0 0 5*a],10,'BaseAngle',60,'Width',1); % plot Z axis arrow
text(-0.01*R0,0,-0.1*a,'0','FontSize',15,'VerticalAlignment','Top','HorizontalAlignment','Right');
text(-0.01*R0,0,5*a,'Z','FontSize',15,'HorizontalAlignment','Right');
arrow([-0.1*R0 0 0],[1.2*(R0+a) 0 0],10,'BaseAngle',60,'Width',1); % plot R axis arrow
text(1.2*(R0+a),0,-0.1*a,'R','FontSize',15,'VerticalAlignment','Top','HorizontalAlignment','Left');

text(R0+a/2,0,a/2,'\theta','FontSize',15,'VerticalAlignment','Middle','HorizontalAlignment','Left');
arrow([R0+0.2*a -0.2*a 0],[R0+0.8*a -0.2*a 0],8,'BaseAngle',60,'Width',1); % plot psi axis arrow
text(R0+0.8*a,-0.9*a,0,'\psi','FontSize',15,'VerticalAlignment','Bottom','HorizontalAlignment','Center');
arrow([R0+1.2*a,0.4*a,0],[R0+1.2*a,-0.4*a,0],8,'BaseAngle',60,'Width',1);
text(R0+1.3*a,0.4*a,0,'\zeta=-\phi','FontSize',15,'VerticalAlignment','Bottom','HorizontalAlignment','Left');

% phi=phimin:pi/20:phimax;
theta=0:pi/20:2*pi;
rmax=a;
phi=0:pi/50:2*pi;
Xf1=[(R0+0*cos(0))*cos(phi),(R0+a*cos(0))*cos(wrev(phi))];
Yf1=[(R0+0*cos(0))*sin(phi),(R0+a*cos(0))*sin(wrev(phi))];
Zf1=0.*Xf1;
Xf2=(R0+rmax*cos(theta))*cos(phimax);
Yf2=(R0+rmax*cos(theta))*sin(phimax);
Zf2=rmax.*sin(theta);
Xf3=(R0+a*cos(pi)).*cos(phi);
Yf3=(R0+a*cos(pi)).*sin(phi);
Zf3=0.*Xf3;
hold on;fill3(Xf1,Yf1,Zf1,'y','facealpha',0.5); % plot S_p
hold on;fill3(Xf2,Yf2,Zf2,'g'); % plot S_t
hold on;fill3(Xf3,Yf3,Zf3,'y','facealpha',0.5); % plot S_p'

hold on; plot3(R0*cos(phimin),R0*sin(phimin),0,'r.','MarkerSize',15);
theta=0:pi/20:2*pi; rmin=0.5*a;
X21=(R0+rmin*cos(theta))*cos(phimin);
Y21=(R0+rmin*cos(theta))*sin(phimin);
Z21=rmin*sin(theta);
hold on; plot3(X21,Y21,Z21,'r','LineWidth',2);
xarwa=R0*cos(phimin);yarwa=R0*sin(phimin);zarwa=rmin;
xarwb=(R0-0.1*a)*cos(phimin);yarwb=(R0-0.1*a)*sin(phimin);zarwb=rmin;
arrow([xarwa,yarwa,zarwa],[xarwb,yarwb,zarwb],10,'BaseAngle',60,'Width',1); % add B_p arrow and text
text(xarwa,yarwa,zarwa,'B_p','FontSize',15,'VerticalAlignment','Middle','HorizontalAlignment','Right');
xp=(R0+a)*cos(phimin-0.3*pi);yp=(R0+a)*sin(phimin-0.3*pi);zp=0;
text(xp,yp,zp,'S_p''','FontSize',15,'VerticalAlignment','Bottom','HorizontalAlignment','Center');

phi=(phimax-0.1*pi):pi/20:(phimax+0.05*pi);
X22=(R0+0*cos(0))*cos(phi);
Y22=(R0+0*cos(0))*sin(phi);
Z22=0.*X22;
hold on; plot3(X22,Y22,Z22,'r','LineWidth',2);
xarwa2=(R0+0*cos(0))*cos(max(phi));yarwa2=(R0+0*cos(0))*sin(max(phi));zarwa2=0;
xarwb2=(R0+0*cos(0))*cos(max(phi)+0.02*pi);yarwb2=(R0+0*cos(0))*sin(max(phi)+0.02*pi);zarwb2=0;
arrow([xarwa2,yarwa2,zarwa2],[xarwb2,yarwb2,zarwb2],10,'BaseAngle',60,'Width',1); % add B_p arrow and text
text(xarwb2,yarwb2,zarwb2,'B_t','FontSize',15,'VerticalAlignment','Bottom','HorizontalAlignment','Left');
xt=(R0+0.8*a*cos(0))*cos(max(phi));yt=(R0+0.8*a*cos(0))*sin(max(phi));zt=0;
text(xt,yt,zt,'S_t','FontSize',15,'VerticalAlignment','Bottom','HorizontalAlignment','Center');
text(0.2*R0,0.2*R0,0,'S_p','FontSize',15,'VerticalAlignment','Bottom','HorizontalAlignment','Center');

% alpha(1);

axis off;
prtstr='poloidal_toroidal_flux_psip_psit';
print('-dpng',[prtstr,'.png']); % output figure to *.png and *.eps
print('-depsc',strcat(prtstr,'.eps'));
% print('-dpdf',strcat(prtstr,'.pdf'));
