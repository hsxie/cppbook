% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-03-13 13:50
% R. B. White's ORBIT code, subroutine spln(ispln,f1,f2,f3,f4,f5,
% f6,f7,f8,f9) in "eqsub.f"
% B-Spline function for f(psi,theta) [B, g, I, ...] and df/dpsi, df/dtheta
% f(x)=f1(j)+f2(j)*h+f3(j)*h^2, h=x-x(j), x(j) <= x < x(j+1)
% uniform x grid
function [f1,f2,f3,f4,f5,f6,f7,f8,f9]=spln2d(ispln,psi,theta,f)
% ispln=1: asymptotic match at pol=0. If ispln=0: original method.
% ccc   perturbation spline calculation, dp = pol - pj, dt = thet - tk
% ccc     grid given by  tk = (k-1)*pi2/lst,  pj = (j-1)*pw/(lsp-1)
% ccc-in domain j,k  f(pol,thet) = f1 + f2*dp + f3*dp**2 + f4*dt + f5*dt*dp
% cccc  + f6*dt*dp**2 + f7*dt**2 + f8*dt**2*dp + f9*dt**2*dp**2
% ccc-  gelg solves bmat(l,m)x(m) = wk2(m), m=1,ndim
% ccc   wk2 contains x after the call
% cccc  bmat stored as vector bmat(lm),lm = l + (m-1)*ndim
% cccc  before calling spln load function f into f1(j,k)=f(pj,tk)
% cccc   after call,  f2(j,k),f3(j,k) etc will be determined

    f1=f';
    f2=0.*f'; f3=0.*f'; f4=0.*f'; f5=0.*f';
    f6=0.*f'; f7=0.*f'; f8=0.*f'; f9=0.*f';

    pw=psi(end)-psi(1);
    lsp=length(psi);
    lst=length(theta);
    lspm = lsp - 1;
    dpx = pw/lspm;
%     dtx = 2*pi/(lst-1);
    dtx = 2*pi/(lst-0);
    
    if(ispln==1)
        % smooth sqrt(f)*df/dp
        % first adjust f1(1,k), use average as on axis quantites.
        f0=0.0;
        for k=1:lst
            f1(1,k)=2.0*f1(2,k)-f1(3,k)-(0.75-sqrt(1.0/32.0))*sqrt(dpx)*f2(1,k);
            f0=f0+f1(1,k);
        end
        f0=f0/lst;

        for k=1:lst
            f1(1,k)=f0;
            % adjust f2(1,k) to enforce continuity of second order derivative
            f2(1,k)=(2.0*f1(2,k)-f1(1,k)-f1(3,k))*8.0/(6.0-sqrt(2.0))/sqrt(dpx);
        end
    else
        for k=1:lst
            f2(1,k) = (10*f1(2,k) - 7*f1(1,k) - 3*f1(3,k))/(4*dpx);
        end
    end
    
    for k = 1:lst
        if(ispln == 1)
            % match f2(2,k) at the first flux surface and smooth the second one
            f2(2,k)=-f2(1,k)*0.5/sqrt(dpx)+(f1(2,k)-f1(1,k))/dpx;
            f1(3,k)=0.4*dpx*f2(2,k)+0.3*f1(4,k)+0.7*f1(2,k);
        end
        for j = (2+ispln):lspm
            jm = j - 1;
            jp = j + 1;
            jpp = min(j + 2,lsp);
            f2(j,k) = -f2(jm,k) + 2*(f1(j,k)-f1(jm,k))/dpx;
            % smooth f1
            if(jp ~= lsp)
                f1(jp,k) = 0.4*dpx*f2(j,k) + 0.3*f1(jpp,k) + 0.7*f1(j,k);
            end
        end
    end
    for k = 1:lst
        f2(lsp,k) = f2(lspm,k);
        % match f3(1,k) at the first flux surface
        if(ispln == 1)
            f3(1,k)=f2(2,k)-f2(1,k)*0.5D0/sqrt(dpx);
        end
        for j = (1+ispln):lspm
            jp = j + 1;
            f3(j,k) = (f2(jp,k)-f2(j,k))/(2*dpx);
        end
        f3(lspm,k) = (f1(lsp,k)-f1(lspm,k)-f2(lspm,k)*dpx)/dpx^2;
    end
% ccc-f1,f2,f3-finished
    
%    ccc   find matrix for f4
    lmax = lst;    
    
    bmat=zeros(lmax,lmax);
    for j=1:lmax
        if(j==lmax)
            jp=1;
        else
            jp=j+1;
        end
        bmat(j,j)=dtx;
%         bmat(j,jp)=dtx;
        bmat(jp,j)=dtx;
    end
    
    
    for j = 1:lspm
        % clear wk2
        wk2=zeros(1,lmax);
        % right hand side
        for k = 1:lst
            kp = k + 1;
            if(k==lst) 
                kp = 1;
            end
            wk2(k) = 2*f1(j,kp) - 2*f1(j,k);
        end

        wk2=wk2/bmat;

        for k = 1:lst
            f4(j,k) = wk2(k);
        end
        for k = 1:lst
            kp = k + 1;
            if(k==lst) 
                kp = 1;
            end
            f7(j,k) = (f4(j,kp) - f4(j,k))/(2*dtx);
        end
    end    
% ccc-f4-f7-finished

% ccc   find matrix for f5
    for j = 1:lspm
        % clear wk2,bmat
        wk2=zeros(1,lmax);

        % ccccc-right hand side
        for k = 1:lst
            kp = k + 1;
            if(k==lst) 
                kp = 1;
            end
            wk2(k) = 2*f2(j,kp) - 2*f2(j,k);
        end

        wk2=wk2/bmat;

        for k = 1:lst
            f5(j,k) = wk2(k);
        end
        for k = 1:lst
            kp = k + 1;
            if(k==lst) 
                kp = 1;
            end
            f8(j,k) = (f5(j,kp) - f5(j,k))/(2*dtx);
        end
    end
% ccc-f5-f8-finished

% ccc   find matrix for f6
    for j = 1:lspm

        % clear wk2,bmat
        wk2=zeros(1,lmax);
        % right hand side
        for k = 1:lst
            kp = k + 1;
            if(k==lst) 
                kp = 1;
            end
            wk2(k) = 2*f3(j,kp) - 2*f3(j,k);
        end


        wk2=wk2/bmat;

        for k = 1:lst
            f6(j,k) = wk2(k);
        end
        for k = 1:lst
            kp = k + 1;
            if(k==lst) 
                kp = 1;
            end
            f9(j,k) = (f6(j,kp) - f6(j,k))/(2*dtx);
        end
    end
% ccc-f6-f9-finished 
    f1=f1'; f2=f2'; f3=f3'; f4=f4'; f5=f5'; f6=f6'; f7=f7'; f8=f8'; f9=f9';
    
end