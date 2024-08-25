% Hua-sheng XIE, 2017-02-23 23:58
% close all; clear; clc;

load phi_th_dr_c.mat;
load phi_th_matrix_c.mat;
load phi_th_ivp_c.mat;
load phi_th_pic_b.mat;
run phi_hd7_dat; 

% phi_pic=phi_pic/phi_pic(floor(length(phi_pic)/2)+1);

% th_dr=tt_dr;

% phi_matrix=smooth(phi_matrix);
% phi_matrix=phi_matrix/phi_matrix(floor(length(phi_matrix)/2));

close all;
figure('unit','normalized','position',[0.01 0.1 0.6 0.5],...
    'DefaultAxesFontSize',15);
subplot(121);
plot(th_hd7,real(phi_hd7),'k-',th_dr,real(phi_dr),':',...
    th_matrix,real(phi_matrix),'--',th_ivp,real(phi_ivp),'-.',...
    th_pic,real(phi_pic),'-','Linewidth',2);
 xlim([-13,20]);
xlabel('\theta');ylabel('Re[\phi]'); ylim([-0.2,1.1]);
subplot(122);
plot(th_hd7,imag(phi_hd7),'k-',th_dr,imag(phi_dr),':',...
    th_matrix,imag(phi_matrix),'--',th_ivp,imag(phi_ivp),'-.',...
    th_pic,imag(phi_pic),'-','Linewidth',2);
xlabel('\theta');ylabel('Im[\phi]'); ylim([-0.2,1.1]); xlim([-13,20]);
hleg=legend('HD7-Dong92, \omega=-0.783+0.335i','mgk-DR, \omega=-0.794+0.336i',...
    'mgk-matrix, \omega=-0.791+0.336i','mgk-ivp, \omega=-0.791+0.335i',...
    'mgk-pic, \omega=-0.785+0.336i',2);legend('boxoff');
set(hleg,'FontSize',10);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 3.0]);
print(gcf,'-dpng','cmp_dong92.png','-r100');
% print(gcf,'-dpdf','cmp_dong92.pdf','-r100');
