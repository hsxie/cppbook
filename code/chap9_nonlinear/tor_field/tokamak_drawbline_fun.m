% Hua-sheng XIE, IFTS-ZJU, huashengxie@gmail.com, 2011-06-26 21:36
function tokamak_drawbline_fun(R,r,n,m)
%     R=5.0;r=1.0;n=2;m=7;
    q=m/n;
    theta=0:pi/1000:m*n*2*pi;
    phi=q.*theta;
    strtitle=['B Field Line of Tokamak, q=m/n=',num2str(m),'/',num2str(n)];
    x=(R+r.*cos(theta)).*cos(phi);
    y=(R+r.*cos(theta)).*sin(phi);
    z=r.*sin(theta);
    plot3(x,y,z);axis equal;xlabel('x');ylabel('y');zlabel('z');
    title(strtitle);
end