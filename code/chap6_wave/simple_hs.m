% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2012-11-16 15:16
% simple example to show half spectral method
% the results match the theory predicts perfectly
function simple_hs
    close all; clear; clc;
    global ua ub k
    ua=1.0; ub=0.2; k=0.2;
    y=3;
    f10=0.1; f20=y*f10;
    
    [T,Y]=ode45(@push,0:0.01:6e2,[f10,f20]);
    tt=T;
    f1=Y(:,1); f2=Y(:,2);
    
    we=[k*(ua+ub), k*(ua-ub)];
    A_w1=f10*abs(y+1)/2; A_w2=f10*abs(y-1)/2; % *f10 ?? need check
    
    h=figure('unit','normalized','Position',[0.01 0.47 0.6 0.45]);
    set(gcf,'DefaultAxesFontSize',15);
    
    subplot(311); plot(tt,real(f1),tt,imag(f1),'LineWidth',2); 
    xlabel('t'); ylabel('f1'); axis tight; grid on;
    title(['k=',num2str(k),', ua=',num2str(ua),', ub=',num2str(ub),...
        ', f1(0)=',num2str(f10),', f2(0)=',num2str(f20)]);   
    subplot(312); plot(tt,real(f2),tt,imag(f2),'LineWidth',2); 
    xlabel('t'); ylabel('f2'); axis tight; grid on;
    title(['theory w_{\pm}=',num2str(we(1)),', ',num2str(we(2)),...
        ', A_{\pm}=',num2str(A_w1),', ',num2str(A_w2)]);
    
    Lt=length(tt); % number of sampling
    dfs=2*pi/(tt(end)-tt(1));
    fs=0:dfs:dfs*(Lt-1);
    f1_ft=fft(real(f1))/Lt*2; % *2 ?? need check
    ifs=30;
    subplot(313); plot(fs(1:ifs),abs(f1_ft(1:ifs)),'LineWidth',2);
    title('simulation frequency');ylabel('Amp');xlabel('\omega');
    xlim([0, fs(ifs)]); grid on;
    Amax=1.5*max(abs(f1_ft(1:ifs))); ylim([0,Amax]);
    hold on; plot([we(1),we(1)],[0,Amax],'r--',[we(2),we(2)],[0,Amax],'r--',...
        'LineWidth',2);

%     fid = fopen ('yt_t.txt','w'); out=[tt';real(f1')];
%     fprintf(fid,'%6.2f%12.8f\n', out);
%     fclose (fid);
end

function dy=push(t,y)
    global ua ub k
    % y --> f1, f2
    dy=zeros(2,1);
    dy(1)=-1i*k*ua*y(1)-1i*k*ub*y(2);
    dy(2)=-1i*k*ub*y(1)-1i*k*ua*y(2);
end
