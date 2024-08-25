% 2013-03-04 11:14
close all; clear; clc;
nNt=4;
Nti=[20,100,100,1000];
figure; set(gcf,'DefaultAxesFontSize',15);
for iNt=1:nNt
    Nt=Nti(iNt); dr=1; x(1)=0;y(1)=0;
    for it=1:Nt
        the=2*pi*rand();
        dx=dr*cos(the);
        dy=dr*sin(the);
        x(it+1)=x(it)+dx;
        y(it+1)=y(it)+dy;
    end
    subplot(2,2,iNt);plot(x,y,'--*',x(1),y(1),'ro',x(end),y(end),'rs','LineWidth',2);
    xymax=1.1*max(max(abs(x),abs(y)));axis equal;
    xlim([-xymax,xymax]);ylim([-xymax,xymax]); 
    xlabel('x'); ylabel('y');
    title(['\Delta{r}=',num2str(dr),', Nt=',num2str(Nt)]);
end

print('-dpng','randwalk2d.png');


