% Hua-sheng XIE, huashengxie@gmail.com, FSC-PKU, 2017-01-12 16:40
% Levy flight rand walk, with Cauthy distribution step size

close all; clear; clc;
nNt=4; icase=3;
Nti=[20,100,100,1000];
figure; set(gcf,'DefaultAxesFontSize',15);
for iNt=1:nNt
    Nt=Nti(iNt); x(1)=0; y(1)=0;
    for it=1:Nt        
        if(icase==1) % fixed step size
            dr=1;
        elseif(icase==2)
            dr=randn(); % Guassian step size
        else
            dr=tan((rand()-0.5)*pi); % Cauthy step size
        end
        the=2*pi*rand();
        dx=dr*cos(the);
        dy=dr*sin(the);
        x(it+1)=x(it)+dx;
        y(it+1)=y(it)+dy;
    end
    subplot(2,2,iNt);
%     plot(x,y,'g--',x,y,'b.',x(1),y(1),'ro',x(end),y(end),'rs','LineWidth',2);
    plot(x,y,'--*',x(1),y(1),'ro',x(end),y(end),'rs','LineWidth',2);
    xymax=1.1*max(max(abs(x),abs(y)));axis equal;
    xlim([-xymax,xymax]);ylim([-xymax,xymax]); 
    xlabel('x'); ylabel('y');
    title(['\Delta{r}=',num2str(dr),', Nt=',num2str(Nt)]);
end

print('-dpng',['randwalk2d_icase=',num2str(icase),'_',num2str(randi(99,1)),'.png']);
