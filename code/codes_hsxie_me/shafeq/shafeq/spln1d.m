% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-03-13 12:30
% R. B. White's ORBIT code, subroutine splin1(f1,f2,f3) in "eqsub.f"
% B-Spline function for f(psi) [T, n, q, ...] and df/dpsi
% f(x)=f1(j)+f2(j)*h+f3(j)*h^2, h=x-x(j), x(j) <= x < x(j+1)
% df/dx=f1(j)+2*h*f3(j)
% uniform x grid
function [f1,f2,f3]=spln1d(x,f)
% comments in "eqsub.f"
% ccc   perturbation spline calculation, dp = pol - pj,
% ccc     grid given by    pj = (j-1)*pw/(lsp-1)
% cccc-in domain j,  f(pol) = f1 + f2*dp + f3*dp**2
% cccc  before calling splin1 load function f into f1(j)=f(pj)
% cccc   after call,  f2(j),f3(j) will be determined
    lsp=length(f);
    lspm = lsp - 1;
    dx=x(2)-x(1);
    f1=f;
    f2=0.*f;
    f3=0.*f;
    % set f2(1) to leave f1(2) unmoved by smoothing
    f2(1) = (10*f1(2) - 7*f1(1) - 3*f1(3))/(4*dx);
    for j=2:lspm
        jm=j-1;
        jp=j+1;
        jpp=min(j+2,lsp);
        f2(j) = -f2(jm) + 2*(f1(j)-f1(jm))/dx;
        % smooth f1
        if(jp~=lsp)
            f1(jp) = 0.4*dx*f2(j)+0.3*f1(jpp) + 0.7*f1(j);
        end
    end
    f2(lsp) = f2(lspm);
    for j = 1:lspm
        jp = j + 1;
        f3(j) = (f2(jp)-f2(j))/(2*dx);
    end
    f3(lspm) = (f1(lsp)-f1(lspm)-f2(lspm)*dx)/dx^2;
end

