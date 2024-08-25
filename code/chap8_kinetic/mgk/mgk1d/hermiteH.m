function a=hermiteH(n,x)
%generate the nth hermite polynomial of x
m=0:floor(n/2);
[q1,q2]=size(m);
s=ndims(x);
[g1,g2]=size(x);
X=repmat(x,[ones(1,s), q2]);
m=permute(repmat(m,[g1,1,g2]),[1,3,2]);
a=factorial(n)*sum((-1).^m./(factorial(m).*factorial(n-2*m)).*(2*X).^(n-2*m),3);
end