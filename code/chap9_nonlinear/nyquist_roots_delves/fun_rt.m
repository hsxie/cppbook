% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-10-07 13:00
% funf.m, give f(z)=D'/D which determined by D(z)
function [N,rt]=fun_rt(za,zb,tol)

if nargin<3, tol=[]; end
if isempty(tol), tol=1e-3; end

% Calculate N, the # of roots in complex domain
intf=fun_sp(0,za,zb,tol);
N=real(round(intf));
if((abs(N-intf)>0.2))
    disp('Accuracy not suffient, need larger nx!!');
    rt=NaN+1i*NaN;
    return;
elseif(N==0)
    disp('N=0 in this domain!!');
    rt=NaN+1i*NaN;
    return;        
end

% Calculate s_p=(1/2i*pi)*\int_c(z^p*f'/f)dz
sp=zeros(N,1); MA=zeros(N,N); Mb=zeros(N,1);
for p=1:N
    sp(p)=fun_sp(p,za,zb,tol);
    Mb(p)=-sp(p);
    MA(p,p)=p;
    for jp=1:N
        ind=jp-p;
        if(ind>=1)
            MA(jp,ind)=sp(p);
        end
    end
end
% sigmap=Mb\MA; % wrong
sigmap=MA\Mb; % calculate simagp, from MA*sigmap=Mb
r0=roots([1;sigmap]);
[ri,jr]=sort(imag(r0),'descend');
rt=r0(jr);
