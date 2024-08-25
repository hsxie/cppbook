% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2017-01-15 02:25

close all; clear; clc;
L=2*pi; N=50; Nt=400;
kk=[0.1,0.2,0.5,1.0,1.5,2.0];
h = figure('Unit','Normalized','position',...
    [0.02 0.1 0.6 0.6],'DefaultAxesFontSize',15);
for jk=1:length(kk)
    K=kk(jk);
    subplot(2,3,jk);
    q=L.*rand(N,1); p=L.*rand(N,1);
    strtitle=['K=',num2str(K)];
    % plot(q,p,'ro','MarkerSize',3);hold on;
    for j=1:Nt
    %	ptmp=p;
        p=p+K.*sin(q);
    %	q=q+ptmp;
        q=q+p;
        q=mod(q+10*L,L); p=mod(p+10*L,L);    
        plot(q,p,'.','MarkerSize',3);hold on;
        xlim([0,L]);ylim([0,L]);
    end
    xlabel('q');ylabel('p');title(strtitle);
end
print(gcf,'-dpng',['standard_mapping.png']);