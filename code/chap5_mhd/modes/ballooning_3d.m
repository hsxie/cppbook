% Hua-sheng XIE, huashengxie@gmail.com, 2015-06-17 12:43
% 3d ballooning mode visulization
close all; clear; clc;

pltc=[0.0  1.0  0.0
    1.0  0.0  0.0
    0.2  0.2  1.0
    0.8 0.8 0.0
    1.0  0.6  0.0
    0.9  0.0  0.9
    0.0  0.8  0.8
    0.0  0.0  0.0
    0.6  0.0  0.0
    0.4  0.7  0.4 
    0.0  0.0  0.5 
    0.6  0.0  0.6 
    0.0  0.5  1.0
    ];

figure('units','normalized','position',[0.02,0.07,0.5,0.8],...
    'DefaultAxesFontSize',15);


% axes('position',[0.08,0.52,0.34,0.34]);

n=10;
r=0:0.002:1.0;
q=1.0+2.*r.^2;
ma=floor(n*min(q))+1;
mb=floor(n*max(q))-1;
m=ma:mb; md=length(m);
% rm=0.*m; qm=0.*m;
qm=m./n;
rm=interp1(q,r,qm);
drc=0.03; rc=0.6; dA=0.1;
% phim=zeros();
for jm=1:md
    phim(jm,:)=exp(-(r-rm(jm)).^2/drc.^2).*exp(-(r-rc).^2/dA.^2);
end

subplot(221);
h=plot(r,q,rm,qm,'rx','Linewidth',2); hold on; ylim([0,4]);
jfs=(1:10)+floor((mb-ma)*0.15);
for jm=jfs
%     plot([rm(jm),rm(jm)],[0,4],'-','LineWidth',2); hold on;
    im=jm-min(jfs)+1;
    hrm(jm)=plot([rm(jm),rm(jm)],[0,4],'--','Color',pltc(im,:),...
        'LineWidth',1); hold on;
    lgd(jm)={['m=',num2str(ma+jm-1)]};
end
text(0.05,3.0,'n=10','Fontsize',13);
xlabel('r/a'); ylabel('q');
% hlgd=legend(hrm(jfs),lgd(jfs),2, 'Fontsize', 8); legend('boxoff');
% %# find handles of lines inside legend that have a non-empty tag
% hLegendLines = findobj(hlgd, 'type', 'line', '-and', '-regexp','Tag','[^'']');
% set(hLegendLines, 'XData', [.45, 0.6]);
% set(hlgd,'position',[0.06,0.75,0.2,0.1]);

% 2.
subplot(222);
% axes('position',[0.54 0.53 0.34 0.34]);
for jm=jfs
    im=jm-min(jfs)+1;
    hphim(jm)=plot(r,phim(jm,:),'-','Color',pltc(im,:),'Linewidth',2); 
    lgd(jm)={['m=',num2str(ma+jm-1)]};
    hold on; %ylim([0,4]);
    plot([rm(jm),rm(jm)],[-0.1,1.2],'--','Color',pltc(im,:),...
        'LineWidth',1); hold on;
end
xlabel('r/a'); ylabel('\delta\phi_m');
hlgd=legend(hphim(jfs),lgd(jfs),2, 'Fontsize', 8); legend('boxoff');
hLegendLines = findobj(hlgd, 'type', 'line', '-and', '-regexp','Tag','[^'']');
set(hLegendLines, 'XData', [.45, 0.6]);
set(hlgd,'position',[0.52,0.73,0.2,0.1]);
axis tight;

% 3. map to 2D
subplot(223);
% axes('position',[0.12 0.09 0.4 0.4]);
t=0:pi/100:2*pi;
[rr,tt]=ndgrid(r,t); % (r,theta)
[col,row]=size(rr);
phi=0.*rr;

for jm=1:md
    m=ma+jm-1;
    phi=phi+repmat(phim(jm,:)',1,row).*exp(-1i*m.*tt);
end

X=1+rr.*cos(tt);
Z=rr.*sin(tt);

% pcolor(X,Z,real(phi)); shading interp; 
contourf(X,Z,real(phi),100,'Linestyle','none');

xlabel('X'); ylabel('Z');

axis equal;
% colormap(hsv); 
% colorbar('location','EastOutside');

%% 4. 3d
% subplot(224);
axes('position',[0.5 0.05 0.45 0.45]);
jr=find(abs(r-rc)==min(abs(r-rc)));
jr=jr(1);
phimc=phim(:,jr);
t3=0.0*pi:0.005*pi:2*pi; p3=-0.5*pi:0.01*pi:1.2*pi;
[tt3,pp3]=meshgrid(t3,p3);
phitp=0.*tt3;
R0=4;
for jm=1:md
    m=ma+jm-1;
    phitp=phitp+phimc(jm).*exp(-1i*(m.*tt3-n.*pp3));
end
xx3=(R0+cos(tt3)).*cos(pp3);
yy3=(R0+cos(tt3)).*sin(pp3);
zz3=sin(tt3);
% contourf(tt3,pp3,real(phitp),100,'Linestyle','none');
% pcolor(tt3,pp3,real(phitp)); shading interp; 
surf(xx3,yy3,zz3,real(phitp),'EdgeColor','none');
axis equal; %mesh off;
box on;

% axis tight;
% axis off;

set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng','ballooning_3d.png');
print(gcf,'-dpdf','-painters','ballooning_3d.pdf');
saveas(gcf,'ballooning_3d.fig','fig');