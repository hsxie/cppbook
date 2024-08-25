% Hua-sheng XIE, 2015-05-10 19:44
close all; clear; clc;
icase=1;
if(icase==1)
    f_fun=@(x,y)-2*pi*pi*sin(pi*x).*sin(pi*y);
    u_fun=@(x,y)sin(pi*x).*sin(pi*y); Lx=2.0; Ly=1.0;
else
    f_fun=@(x,y)-2*((1-6*x.^2).*y.^2.*(1-y.^2)+(1-6*y.^2).*x.^2.*(1-x.^2));
    u_fun=@(x,y)-(x.^2-x.^4).*(y.^2-y.^4); Lx=1.0; Ly=1.0;
end
eps=1e-8; nx=128; ny=64; dx=Lx/(nx+0); dy=Ly/(ny+0); w=1.0; 
dx2r=1.0/(dx*dx); dy2r=1.0/(dy*dy);
[xx,yy]=ndgrid(0:dx:Lx,0:dy:Ly);
ff=f_fun(xx,yy); uu=0.*xx; ncout=0; maxdu=2*eps;
while (maxdu>eps && ncout<10000)
    tmpu=uu;
    for i=2:nx
        for j=2:ny 
            uu(i,j)=(1.0-w)*uu(i,j)-w/(2.0*dx2r+2.0*dy2r)*(ff(i,j)...
                -uu(i+1,j)*dx2r-uu(i-1,j)*dx2r...
                -uu(i,j+1)*dy2r-uu(i,j-1)*dy2r);
        end
    end
    maxdu=max(max(abs(tmpu-uu))); ncout=ncout+1;
end
figure('DefaultAxesFontSize',15);
subplot(221); mesh(xx,yy,ff); xlabel('x');  axis tight;
ylabel('y'); zlabel('f'); box on;
subplot(222); mesh(xx,yy,uu); xlabel('x');  axis tight;
ylabel('y'); zlabel('u_{numerical}'); box on;
subplot(223); mesh(xx,yy,u_fun(xx,yy));  axis tight;
xlabel('x'); ylabel('y'); zlabel('u_{exact}'); box on;
subplot(224); mesh(xx,yy,uu-u_fun(xx,yy));  axis tight;
xlabel('x'); ylabel('y'); zlabel('u_{numerical}-u_{exact}'); box on;

