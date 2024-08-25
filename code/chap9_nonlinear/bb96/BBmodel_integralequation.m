% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2011-08-03 11:42,
% modified 2011-08-20 16:31, 2011-11-20 18:47
% BBmodel_heeter.m, 2nd version.
% Here we use the method of Heeter1999, from Lilley2009 PhD thesis P131, 
% which can reduce the computational intensity of the simulation from
% O(nt^3) to O(nt^2). 
% Typical run time: nt=3000, t=9s; nt=12000, t=145s.
% method=1, change two-dimension-array S(j,k) to two one-dimension-array 
% S(k) and S1(k) in case of large nt (e.g. 2*10^4); method != 1, keep 
% S(j,k) be two dimensions, this is straightfoward and can be faster by a
% slight, but may out of memory when nt is too large.
close all;clear;clc;
method=1;

h=figure('unit','normalized','position',[0.02,0.1,0.6,0.7],...
    'DefaultAxesFontSize',12);
for jplt=1:4
    
runtime=cputime;

if(jplt==1)
    % a is alpha, b is beta, v is nu
    nt=3000;dt=0.01; A0=0.1; a=0;b=0.0;v=4.31;
elseif(jplt==2)
    nt=12000;dt=0.01; A0=0.1; a=0;b=0.0;v=2.18;
elseif(jplt==3)
    nt=20000;dt=0.005; A0=0.1; a=0;b=0.0;v=1.28;
else
    nt=1000;dt=0.01; A0=0.1; a=0;b=0.0;v=1.15;
end

c1=-dt^5/2.0;c2=b*dt;c3=(v*dt)^3;c4=(a*dt)^2;
A=repmat(A0,1,nt+1);
tt=repmat(0.0,1,nt+1);
switch method
    case 1
        S=repmat(0.0,1,floor(nt/2)+1);
        S1=repmat(0.0,1,floor(nt/2)+1); % a temp array
        for j=0:nt-1
            IntAS=0.0;
            kmax=floor(j/2);
            for k=1:kmax
                if(k==kmax)
                    S(k)=0.0; % S(j,k)
                    for l=0:(j-2*k)
                        S(k)=S(k)+exp(c2*(l-j)+c3*k^2*(4.0*k/3.0+l-j)+1i*c4*k*(j-l-k))*A(l+k+1)*A(l+1);
                        S1(k)=S(k);
                    end
                else
                    S(k)=exp(-c2-c3*k^2+1i*c4*k)*S1(k)+exp(-2*k*c2-2*c3*k^3/3+1i*c4*k^2)*A(j-k+1)*A(j-2*k+1);
                    S1(k)=S(k);
                end
                IntAS=IntAS+k^2*A(j-k+1)*S(k); % Note: k^2 is missed in Lilley2009 (6.18).
            end
            A(j+1+1)=A(j+1)+A(j+1)*dt+c1*IntAS;
            tt(j+1+1)=j*dt;
        end
    otherwise
        S=repmat(0.0,nt+1,round(nt/2)+1);
        for j=0:nt-1
            IntAS=0.0;
            kmax=round(j/2);
            for k=1:kmax
                if(k==kmax)
                    S(j+1,k)=0.0;
                    for l=0:(j-2*k)
                        S(j+1,k)=S(j+1,k)+exp(c2*(l-j)+c3*k^2*(4.0*k/3.0+l-j)+1i*c4*k*(j-l-k))*A(l+k+1)*A(l+1);
                    end
                else
                    S(j+1,k)=exp(-c2-c3*k^2+1i*c4*k)*S(j-1+1,k)+exp(-2.0*k*c2-2.0*c3*k^3/3.0+1i*c4*k^2)*A(j-k+1)*A(j-2*k+1);
                end
                IntAS=IntAS+k^2*A(j-k+1)*S(j+1,k); % Note: k^2 is missed in Lilley2009 (6.18).
            end
            A(j+1+1)=A(j+1)+A(j+1)*dt+c1*IntAS;
            tt(j+1+1)=j*dt;
        end
end
runtime=cputime-runtime;
% plot(tt,abs(A));ylabel('|A|');
subplot(2,2,jplt);
plot(tt,A,'linewidth',2);ylabel('A');
xlabel('\tau');xlim([min(tt),max(tt)]);
if(jplt==4)
    ylim([-100,1e2]);
end
title(['Run time ',num2str(runtime,4), 's',10,...
    'A(0) =',num2str(A0),', \Delta t =',num2str(dt),', \alpha =',...
    num2str(a),', \beta =',num2str(b),', \nu =',num2str(v)]);
end
