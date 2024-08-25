% Hua-sheng XIE, 2014-07-12 23:02
% test fft to minus/plus directions

close all; 
clear; clc;

dt=0.05; dx=0.02;
t=0:dt:50;
x=0:dx:10;

k=2; w1=3.2; w2=-5.1; A1=5; A2=3; %%
[tt,xx]=meshgrid(t,x);
yy=A1*cos(-k*xx-w1*tt)+A2*cos(2*k*xx-w2*tt);
y=A1*cos(w1*t)+A2*cos(w2*t);

nt=length(t); nx=length(x);

yf=fftshift(fft(y));
yyf=fftshift(fft2(yy));
tf0=2*pi/dt.*linspace(-0.5,0.5,nt);
tf=tf0;
dtf=tf0(2)-tf0(1);


h = figure('Unit','Normalized','position',...
    [0.02 0.2 0.5 0.6],'DefaultAxesFontSize',15);

subplot(221);plot(t,y,'LineWidth',2);
xlim([0,10]); xlabel('t');
title(['(a) y(t), k=',num2str(k),', \omega_1=',num2str(w1),...
    ', \omega_2=',num2str(w2)]);

subplot(222);
plot(tf,real(yf),tf,imag(yf),'--','LineWidth',2);
xlim([-10,10]);
xlabel('\omega');
% axis tight;
title('(b) y(t) power spectral');

%%
subplot(223);
pcolor(xx,tt,yy); shading interp;
% surf(xx,tt,yy);
xlim([0,5]); ylim([0,10]);
xlabel('x'); ylabel('t');
title(['(c) y(x,t), A_1=',num2str(A1),', A_2=',num2str(A2)]);

subplot(224);
% tf=1/dt.*linspace(0,1,nt);
kf=2*pi/dx.*linspace(-0.5,0.5,nx);
[ttf,kkf]=meshgrid(tf,kf);
% surf(kkf,ttf,abs(yyf));
pcolor(kkf,ttf,abs(yyf));
% pcolor(kkf,ttf,real(yyf));
% pcolor(kkf,ttf,imag(yyf));
shading interp;
xlim([-5,5]); ylim([-10,10]);
title('(d) \omega v.s. k'); 
xlabel('k'); ylabel('\omega');

print(gcf,'-dpng',['tst_fft2_wk_k=',num2str(k),'_w1=',num2str(w1),...
    '_w2=',num2str(w2),'_A1=',num2str(A1),',A2=',num2str(A2),'.png']);

