[XX,YY]=meshgrid(-1:0.01:1);
[Q,R]=cart2pol(XX,YY);
R(find(R>1))=NaN;
h=figure('unit','normalized','position',[0.02,0.1,0.6,0.45],...
    'DefaultAxesFontSize',12);
for m=0:7
    subplot(2,4,m+1);
    pcolor(XX,YY,exp(-50.*(R-0.5).^2).*cos(m*Q));
    title(['m=',num2str(m)]);
    shading interp; box on; axis equal;
end