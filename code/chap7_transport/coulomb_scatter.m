% Hua-sheng XIE, 2017-01-05 20:28
close all; clear; clc;

xv0=[-2.0, 0.5, 1.0, 0.0;
    -2.0, 1.0, 1.0, 0.0;
    -2.0, 1.5, 1.0, 0.0;
    -2.0, 2.0, 1.0, 0.0;
    -2.0, 1.0, 0.5, 0.0;
    -2.0, 1.0, 1.5, 0.0;
    -2.0, 1.0, 2.0, 0.0;
    -3.0, 1.5, 1.0, 0.0];

h=figure('unit','normalized','Position',[0.01 0.17 0.4 0.45],...
    'DefaultAxesFontSize',15);

cmap=colormap('jet');
ncmap=length(cmap);
clra=min(xv0(:,3)); clrb=max(xv0(:,3));

for jp=1:length(xv0)
%     x=-2.0; y=0.5; vx=1.0; vy=0.0;
    x=xv0(jp,1);y=xv0(jp,2);vx=xv0(jp,3);vy=xv0(jp,4);
    nt=1200; dt=0.005;
    tmp=4;
    dt=dt/tmp; nt=nt*tmp;
    xx=[];yy=[];vxx=[];vyy=[];tt=[];

    plot([x,-x],[0,0],'m--',[x,-x],[y,y],'k:',0,0,'ro',...
        'MarkerSize',5,'MarkerFaceColor','r'); hold on;
    for it=1:nt
        r3=(sqrt(x^2+y^2))^3;
        Fx=-x/r3;
        Fy=-y/r3;
        x=x+vx*dt;
        y=y+vy*dt;
        vx=vx+Fx*dt;
        vy=vy+Fy*dt;
        tt=[tt,it*dt];
        xx=[xx,x];
        yy=[yy,y];
    end
    plot(xx,yy,'linewidth',2,'color',cmap(floor((xv0(jp,3)-clra)/(clrb-clra)*(ncmap-1))+1,:)); hold on;
end
%%
dclr=(clrb-clra)/7;
colorbar('YTickLabel',{num2str(clra+1*dclr,2),num2str(clra+2*dclr,2),...
    num2str(clra+3*dclr,2),num2str(clra+4*dclr,2),num2str(clra+5*dclr,2),...
    num2str(clra+6*dclr,2)});
xlabel('x');ylabel('y');title('coulomb scatter');axis tight;
print(gcf,'-dpng','coulomb_scatter.png');
