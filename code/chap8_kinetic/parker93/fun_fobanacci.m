% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-08-06 13:23
% Uniform loading in x and with non-random loading of Gaussian in v with
% Fobanacci numbers 
% 16-10-24 10:37
function [x,v,N]=fun_fobanacci(M)
% M=23;
alpha=ones(1,M+1);
for m=3:M+1
    alpha(m)=alpha(m-1)+alpha(m-2);
end
alpha(1)=[];
N=alpha(M);
vth=1.0;
x=zeros(1,N); y=0.*x; v=0.*x;
method=2; %
for i=1:N
    x(i)=(2*i-1)/(2*alpha(M));
    y(i)=alpha(M-1)*x(i);
    if(method==1)
        v(i)=sqrt(2)*vth*erfcinv(2*x(i));
    else
        y(i)=mod(y(i),1);
        v(i)=sqrt(2)*vth*erfcinv(2*y(i));
    end
end
v(v==Inf)=0.0;

