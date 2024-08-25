% Hua-sheng XIE, 2016-10-02 15:03
% Ref: Cummings1994 thesis, sec 4.1 for 1D gyroaverage weight and point
close all; clear; clc;
nring=8*2;

figure('Unit','Normalized','position',...
    [0.02 0.1 0.3 0.5],'DefaultAxesFontSize',15);
t=(0:pi/50:2*pi)-pi/1;
x=cos(t); y=sin(t);
plot(x,y,'g--','Linewidth',2); hold on;
axis equal;
xlabel('x'); ylabel('y');


for j=1:nring
    t2d=(2*j-1)/nring*pi;
    x2d=cos(t2d); y2d=sin(t2d);
    plot(x2d,y2d,'kx','Linewidth',2,'Markersize',8); hold on;
    
    t2d=(2*j)/nring*pi;
    x2d=cos(t2d); y2d=sin(t2d);
    plot(x2d,y2d,'ro','Linewidth',2,'Markersize',5); hold on;
end

wgt=[];xxt=[];
for j=1:nring
    t1=(2*j-1)/nring*pi; t2=(2*j+1)/nring*pi;
    xt=(sin(t2)-sin(t1))/(t2-t1);
    wgt=[wgt,1/nring];
    xxt=[xxt,xt];
%     y1=cos(t1); y2=cos(t1);
    plot(xt,-1,'b^','Linewidth',2,'Markersize',5);
    hold on;
end
title(['nring=',num2str(nring)]);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 4]);
print(gcf,'-dpng',['plt_wgt_nring=',num2str(nring),'.png'],'-r100');

