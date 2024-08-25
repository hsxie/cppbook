% plot Takeda1991 Fig. 2.5
% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2013-01-10 13:54
close all; clear; clc;
R0=1.0; a=0.3; delta=0.5; kappa=1.2; % parameters

theta=0:pi/20:2*pi; % cal Tokamak shape
r=R0+a.*cos(theta+delta.*sin(theta));
z=kappa*a.*sin(theta);
ymax=kappa*a; xmax=R0+a;

sftind=find(z==ymax);
xsft=r(sftind(1)); % find x_shift position

figure; set(gcf,'DefaultAxesFontSize',15);
plot(r,z,'r','LineWidth',2); axis equal;hold on; % plot Tokamak shape

xlim([-0.2*xmax,1.4*xmax]);ylim([-1.5*ymax,1.5*ymax]);
arrow([-0.1*xmax 0],[1.2*xmax 0],10,'BaseAngle',60,'Width',1); % plot x axis arrow
text(1.15*xmax,-0.01*ymax,'R','FontSize',15,'VerticalAlignment','Top');
arrow([0 -1.2*ymax],[0 1.5*ymax],10,'BaseAngle',60,'Width',1); % plot y axis arrow
text(-0.01*R0,1.3*ymax,'Z','FontSize',15,'HorizontalAlignment','Right');
text(-0.01*R0,-0.01*ymax,'0','FontSize',15,'HorizontalAlignment',...
    'Right','VerticalAlignment','Top'); % text (0,0) point
axis off;

text(R0,0,'R_0','FontSize',15,'VerticalAlignment','Top'); % plot x=R0 line
hold on;plot([R0,R0],[0,1.3*ymax],'--','LineWidth',2);
hold on;plot(R0,0,'g.','MarkerSize',15);

text(-0.01*R0,ymax,'\kappa{a}','FontSize',15,...
    'HorizontalAlignment','Right'); % plot y=kappa*a line
hold on;plot([0,xmax],[ymax,ymax],'--','LineWidth',2);

text((xsft+R0)/2,1.25*ymax,'\Delta','FontSize',15,...
    'HorizontalAlignment','Center'); % plot x=x_shift line
text(xsft,1.25*ymax,'\rightarrow','FontSize',15,...
    'HorizontalAlignment','right');
text(R0,1.25*ymax,'\leftarrow','FontSize',15,...
    'HorizontalAlignment','left');
hold on;plot([xsft,xsft],[0,1.3*ymax],'--','LineWidth',2);

text(R0+a/2,1.15*ymax,'a','FontSize',15,...
    'HorizontalAlignment','Center'); % plot x=xmax line
text(R0,1.15*ymax,'\leftarrow','FontSize',15,...
    'HorizontalAlignment','Left');
text(R0+a,1.15*ymax,'\rightarrow','FontSize',15,...
    'HorizontalAlignment','Right');
hold on;plot([xmax,xmax],[0,1.3*ymax],'--','LineWidth',2);

text(0.01*R0,0.8*ymax,['R=R_0+acos(\theta+\delta sin\theta)',10,...
    'Z=\kappa asin\theta'],'FontSize',15,...
    'VerticalAlignment','Top'); % legend equation
text(0.01*R0,-0.8*ymax,['Ellipticity: \kappa',10,...
    'Triangularity: \Delta=asin\delta'],'FontSize',15); % legend parameter

tstr=['Here: R_0=',num2str(R0),', a=',num2str(a),...
    ', \kappa=',num2str(kappa),', \delta=',num2str(delta)];
title(tstr);

prtstr=['Tokamak_Shape,R0=',num2str(R0),',a=',num2str(a),...
    ',kappa=',num2str(kappa),',delta=',num2str(delta)];
print('-dpng',[prtstr,'.png']); % output figure to *.png and *.eps
print('-depsc',strcat(prtstr,'.eps'));
