% Hua-sheng XIE, IFTS-ZJU, 2011-07-02, 15:36
% standardmap_fun.m, p(j+1)=p(j)+Ksin(q(j)), q(j+1)=q(j)+p(j)
function standardmap_fun(K,Npoints,Nstep,plotpq0,pperiod,initial,hAxes1)
%     close all;clear;clc;
%     K=4.0;N=200;maxj=500;
    set(gcf,'currentaxes',hAxes1);
%     set(hAxes1,'color',[0.0 0.0 0.0]);
    L=2*pi;N=Npoints;
    if(initial==1)
        q=L.*rand(N,1);
        p=L.*rand(N,1);
    elseif(initial==2)
        p=linspace(0,L,N);
        q=0.*p;
    else
        q=linspace(0,L,N);
        p=0.*q;
    end
    strtitle=['k=',num2str(K)];
%     figure;
    if(plotpq0)
        plot(q,p,'ro','MarkerSize',3);hold on;
    end
    xlabel('q');ylabel('p');title(strtitle);
%     scatter(q,p,3);
    xlim([0,L]);ylim([0,L]);hold on;
    for j=1:Nstep
%         ptmp=p;
        p=p+K.*sin(q);
    %     q=q+ptmp; 
        q=q+p;
        % bc
        q=q./L+100.0;
        q=L.*(q-floor(q));
        if(pperiod)
            p=p./L+100.0;
            p=L.*(p-floor(p));
        end
%         set(hAxes1,'color',[0.0 0.0 0.0]);
        hold on;
        plot(q,p,'.','MarkerSize',3);
%         scatter(q,p,1,'.');
        pause(0.00001);
    end
end