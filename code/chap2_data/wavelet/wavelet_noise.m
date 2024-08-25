close all;clear;clc;
t =0:0.01:10;
yt0=sin(pi.*t)+0.5.*cos(2.0*pi.*t);
ytn=yt0+2.0.*(rand (1,length(t))-0.5);
% [XD,CXD,LXD] = wden(X,TPTR,SORH,SCAL,N,'wname');
% wname: 'db1' or 'haar', 'db2', ...'sym2', ... , 'sym8', ... 
ytd=wden(ytn,'minimaxi','s','one',5,'db3');
figure; set (gcf,'DefaultAxesFontSize',12);
subplot(211);plot(t,ytn,'linewidth',2);title('noise data');ylim([-2.5,2.5]);
subplot(212);plot(t,yt0,t,ytd,'r','linewidth',2);ylim([-2.5,2.5]);
legend('original data','de-noising data');legend('boxoff');

% fid = fopen ('yt_t .txt ', 'w '); out =[t;yt ];
% fprintf (fid , '%6.2 f %12.8 f\n', out);
% fclose (fid);

% N=1000; t=1:N; x=sin(0.05*t);
% load noissin;
% ns=noissin;
% subplot(311);plot(t,x);title('original data');
% subplot(312);plot(ns);title('noise data');
% xd=wden(ns,'minimaxi','s','one',5,'db3');
% subplot(313);plot(t,xd);title('de-noising data');
