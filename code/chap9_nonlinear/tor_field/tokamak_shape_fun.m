% Hua-sheng XIE @ IFTS-ZJU, huashengxie@gmail.com, 2011-04-19 11:01
% R.B.White, The Theory of Toraidally Confined Plasmas, 2001
% P143, (4.169)-(4.170)
% x0 is the major radius, a the minor radius, kappa the ellipticity, delta
% the triangularity, and b the indentation, which yields a bean shape.
function tokamak_shape_fun(x0,kappa,delta,b0)
%    clear;clc;close all;
%    x0=5.0;kappa=1.4;delta=0.0;
    hold off;
    N=200;theta=1e-4*x0:pi/N:2*pi;[m,n]=size(theta);
    xtmp=x0.*theta./theta;ztmp=0.*theta;
    for a=0.1:0.1:1.0
        b=b0*a;
        str=['\delta =',num2str(delta),10,'\kappa =',num2str(kappa),10,...
            'b =',num2str(b/a),'a'];
        x=x0-b+(a+b.*cos(theta)).*cos(theta+delta.*sin(theta));
        z=kappa.*a.*sin(theta);
        plot(x,z);axis equal;hold on;
        for j=1:n
            if(mod(j,floor(n/20))==0)
                plot([x(j),xtmp(j)],[z(j),ztmp(j)],'-');hold on;
            end
        end
        xtmp=x;ztmp=z;
    end
    text(x0+0.7*a,0.95*a,str);title('Tokamak Cross Section');
    xlabel('R');ylabel('Z');
end