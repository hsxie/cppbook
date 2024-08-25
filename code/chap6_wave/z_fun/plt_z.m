% Hua-sheng XIE, 2017-04-06 23:02
% Plot dispersion function with and without analytical continuation
close all; clear; clc;
zeta=@(x)faddeeva(x)*1i*sqrt(pi);
zeta_up=@(x)Z_fun(x);

fz=@(z)zeta(z);
fzp=@(z)-2*(1+z.*zeta(z));
fz_up=@(z)zeta_up(z);

%% 1. 
h=figure('unit','normalized','position',[0.02,0.1,0.55,0.5],...
    'DefaultAxesFontSize',12);
x=-6:0.01:6; y=-3:0.01:3;
subplot(221); z=x;
plot(x,real(fz(z)),'-',x,imag(fz(z)),'--',x,real(fzp(z)),'-.',...
    x,imag(fzp(z)),':','linewidth',2);
% legend('Re(Z)','Im(Z)','Re(Zp)','Im(Zp)',1); legend('boxoff'); 
grid on; xlabel('x (y=0)'); axis tight;
subplot(222); z=1i*y;
plot(y,real(fz(z)),'-',y,imag(fz(z)),'--',y,real(fzp(z)),'-.',...
    y,imag(fzp(z)),':','linewidth',2);
legend('Re(Z)','Im(Z)','Re(Zp)','Im(Zp)',1); legend('boxoff'); 
grid on; xlabel('y (x=0)'); axis tight;
subplot(223); z=x-3i;
plot(x,real(fz(z)),'-',x,imag(fz(z)),'--',x,real(fzp(z)),'-.',...
    x,imag(fzp(z)),':','linewidth',2);
% legend('Re(Z)','Im(Z)','Re(Zp)','Im(Zp)',1); legend('boxoff'); 
grid on; xlabel('x (y=-3)'); axis tight;
subplot(224); z=1i*y+6;
plot(y,real(fz(z)),'-',y,imag(fz(z)),'--',y,real(fzp(z)),'-.',...
    y,imag(fzp(z)),':','linewidth',2);
% legend('Re(Z)','Im(Z)','Re(Zp)','Im(Zp)',1); legend('boxoff'); 
grid on; xlabel('y (x=6)'); axis tight;


%% 2. 
% close all;
h=figure('unit','normalized','position',[0.02,0.1,0.45,0.6],...
    'DefaultAxesFontSize',12);

x=-10:0.1:10; y=-10:0.1:10;
[xx,yy]=meshgrid(x,y);
zz=xx+1i*yy;
subplot(221);
contourf(xx,yy,real(fz(zz))); xlabel('x'); ylabel('y');
title('Re(Z), contourf'); axis equal; axis tight;
subplot(222);
contourf(xx,yy,imag(fz(zz))); xlabel('x'); ylabel('y');
title('Im(Z), contourf'); axis equal; axis tight;
subplot(223);
imagesc(x,y,real(fz(zz))); xlabel('x'); ylabel('y'); caxis([-1 1]);
title('Re(Z), imagesc');axis xy square; 
subplot(224);
imagesc(x,y,imag(fz(zz))); xlabel('x'); ylabel('y'); caxis([-1 1]);
title('Re(Z), imagesc');axis xy square; 

%% 3. 
% close all;
h=figure('unit','normalized','position',[0.02,0.1,0.6,0.6],...
    'DefaultAxesFontSize',12);

x=-2.02:0.1:2; y=-2.02:0.1:2;
[xx,yy]=meshgrid(x,y);
zz=xx+1i*yy;
subplot(221);surfc(xx,yy,real(fz(zz)));axis tight; colorbar;
xlabel('x'); ylabel('y'); title('(a) Re(Z), w/ analytical continuation');
subplot(222);surfc(xx,yy,imag(fz(zz)));axis tight; colorbar;
xlabel('x'); ylabel('y'); title('(b) Im(Z), w/ analytical continuation');
subplot(223);surfc(xx,yy,real(fz_up(zz)));axis tight; colorbar;
xlabel('x'); ylabel('y'); title('(c) Re(Z), w/o analytical continuation');
subplot(224);surfc(xx,yy,imag(fz_up(zz)));axis tight; colorbar;
xlabel('x'); ylabel('y'); title('(d) Im(Z), w/o analytical continuation');
