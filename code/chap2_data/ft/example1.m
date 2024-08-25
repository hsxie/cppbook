t=0:0.01:10;
yt=0.1.*sin(15.*2.*pi.*t+0.1*pi).*exp(0.1.*t)-0.4.*cos(5.*2.*pi.*t...
    +0.7*pi)+0.8.*sin(2.*2.*pi.*t-0.1*pi)+rand(1,length(t));
figure; set(gcf,'DefaultAxesFontSize',15);
plot(t,yt,'linewidth',2);title('yt-t');xlabel('t/s');ylabel('yt');
fid = fopen('yt_t.txt', 'w');
out=[t;yt];
fprintf(fid, '%6.2f %12.8f\n', out);
fclose(fid);