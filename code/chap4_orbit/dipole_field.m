% Hua-sheng XIE, huashengxie@gmail.com, 2015-06-12 13:41
% Plot dipole field lines
% 2017-01-03 13:31 Add B=B_dip+ B_T\hat{x}
close all; clear; clc;

B0=1.0; Re=1.0;
BT=0.01*B0; deltaz=0.1*Re;

fr=@(x,y,z)sqrt(x.^2+y.^2+z.^2);
fBx=@(x,y,z)-B0*Re^3*3*x.*z./fr(x,y,z).^5+BT*tanh(-z/deltaz);
fBy=@(x,y,z)-B0*Re^3*3*y.*z./fr(x,y,z).^5;
fBz=@(x,y,z)-B0*Re^3*(2*z.^2-x.^2-y.^2)./fr(x,y,z).^5;
fB=@(x,y,z)sqrt(fBx(x,y,z).^2+fBy(x,y,z).^2+fBz(x,y,z).^2);

figure('units','normalized','position',[0.02,0.1,0.7,0.5],...
    'DefaultAxesFontSize',15);

% ri=[0.1, 0.0, 0.0;
%     0.3, 0.0, 0.0;
%     0.7, 0.0, 0.0;
%     1.3, 0.0, 0.0;
%     1.9, 0.0, 0.0;
%     2.8, 0.0, 0.0;
%     4.8, 0.0, 0.0;
%     6.8, 0.0, 0.0;
%     0.1, 0.0, 3.0;
%     0.5, 0.0, 3.0;
%     1.0, 0.0, 3.0;
%     1.8, 0.0, 3.0;
%     0.1, 0.0, -3.0;
%     0.5, 0.0, -3.0;
%     1.0, 0.0, -3.0;
%     1.8, 0.0, -3.0;
%     ];
ri=[0.1, 0.0, 0.0;
    0.3, 0.0, 0.0;
    0.7, 0.0, 0.0;
    1.3, 0.0, 0.0;
    1.9, 0.0, 0.0;
    2.8, 0.0, 0.0;
    4.8, 0.0, 0.0;
    6.8, 0.0, 0.0;
    0.1, 0.0, 2.0;
    0.5, 0.0, 2.0;
    1.0, 0.0, 2.0;
    1.8, 0.0, 2.0;
    0.1, 0.0, -2.0;
    0.5, 0.0, -2.0;
    1.0, 0.0, -2.0;
    1.8, 0.0, -2.0;
    ];
ri=[ri;-ri];

ri=ri.*Re;

xi=ri(:,1); yi=ri(:,2); zi=ri(:,3);

npl=length(xi);
for theta=0:30:150;
    
    for ipl=1:npl
    %     rnd=rand();
    %     x0=5*Re*rnd; y0=5*Re*rnd; z0=0.0;

%         x0=xi(ipl); y0=yi(ipl); z0=zi(ipl);
        x0=xi(ipl)*cos(theta*pi/180); y0=xi(ipl)*sin(theta*pi/180); z0=zi(ipl);
        dl=0.02; x=[]; y=[]; z=[];
        x(1)=x0; y(1)=y0; z(1)=z0;
    %     subplot(121); plot3(x0,y0,z0,'r+'); hold on;
    %     subplot(122); plot(x0,z0,'r+'); hold on;
        for pm=[-1,1]
            for it=1:600
                if(fr(x(it),y(it),z(it))<0.01*Re || abs(z(it))>8*Re || x(it)<-5.0*Re || x(it)>10*Re)
                    x(it)=NaN; y(it)=NaN; z(it)=NaN;
                end
                x(it+1)=x(it)+dl*fBx(x(it),y(it),z(it))/fB(x(it),y(it),z(it))*pm;
                y(it+1)=y(it)+dl*fBy(x(it),y(it),z(it))/fB(x(it),y(it),z(it))*pm;
                z(it+1)=z(it)+dl*fBz(x(it),y(it),z(it))/fB(x(it),y(it),z(it))*pm;
            end
            subplot(121); plot3(x,y,z,'Linewidth',2); hold on; box on;
            if(theta==0)
                subplot(122); plot(x,z,'Linewidth',2); hold on; box on;
            end
        end
    end
end
%%
subplot(121); xlabel('x'); ylabel('y'); zlabel('z'); axis equal; axis tight;
subplot(122); xlabel('x'); ylabel('z'); axis equal; axis tight;
title(['B_T/B_0=',num2str(BT/B0),', \delta/R_e=',num2str(deltaz/Re)]);

set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng','dipole_field.png');
print(gcf,'-dpdf','-painters','dipole_field.pdf');
saveas(gcf,'dipole_field.fig','fig');
