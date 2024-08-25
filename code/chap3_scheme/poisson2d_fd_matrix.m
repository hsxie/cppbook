% Hua-sheng XIE, 2015-05-11 16:23
close all; clear; clc;
icase=1;
if(icase==1)
    f_fun=@(x,y)-2*pi*pi*sin(pi*x).*sin(pi*y);
    u_fun=@(x,y)sin(pi*x).*sin(pi*y); Lx=2.0; Ly=1.0;
else
    f_fun=@(x,y)-2*((1-6*x.^2).*y.^2.*(1-y.^2)+(1-6*y.^2).*x.^2.*(1-x.^2));
    u_fun=@(x,y)-(x.^2-x.^4).*(y.^2-y.^4); Lx=1.0; Ly=1.0;
end
nx=8*128/4; ny=4*64/2; dx=Lx/(nx+1); dy=Ly/(ny+1);
ha=1.0/(dx*dx); hb=1.0/(dy*dy); hc=ha+hb;
[xx,yy]=ndgrid(0:dx:Lx,0:dy:Ly);
ff=f_fun(xx,yy); uu=0.*xx;

rhs=ff(1:nx,1:ny);
II=1:nx; JJ=1:ny;

% adjust the rhs to include boundary terms:
% rhs(:,1) = rhs(:,1) - usoln(Iint,1)/h^2;
% rhs(:,ny) = rhs(:,ny) - usoln(Iint,ny+2)/h^2;
% rhs(1,:) = rhs(1,:) - usoln(1,Jint)/h^2;
% rhs(nx,:) = rhs(nx,:) - usoln(nx+2,Jint)/h^2;

% convert the 2d grid function rhs into a column vector for rhs of system:
F = reshape(rhs,nx*ny,1);

% form matrix A:
I = speye(ny); I2 = speye(nx);
e = ones(nx,1);
D = spdiags([ha*e -2*hc*e ha*e],[-1 0 1],nx,nx);
S = spdiags([hb*e hb*e],[-1 1],ny,ny);
A = (kron(I,D) + kron(S,I2));

if(nx*ny<1000)
 A_ful=full(A); I_ful=full(I); e_ful=full(e); S_ful=full(S); D_ful=full(D);
end
% Solve the linear system:
uvec = A\F;  
uu(II,JJ) = reshape(uvec,nx,ny);


figure('DefaultAxesFontSize',15);
subplot(221); mesh(xx,yy,ff); xlabel('x');  axis tight;
ylabel('y'); zlabel('f'); box on;
subplot(222); mesh(xx,yy,uu); xlabel('x');  axis tight;
ylabel('y'); zlabel('u_{numerical}'); box on;
subplot(223); mesh(xx,yy,u_fun(xx,yy));  axis tight;
xlabel('x'); ylabel('y'); zlabel('u_{exact}'); box on;
subplot(224); mesh(xx,yy,uu-u_fun(xx,yy));  axis tight;
xlabel('x'); ylabel('y'); zlabel('u_{numerical}-u_{exact}'); box on;

