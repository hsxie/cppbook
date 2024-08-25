% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2016-10-08 13:24
% Divide the original domain to subdomains, with which contain at most M
% zeros, e.g., M=1
function [N,domain]=fun_divide_domain(za,zb,M,tol)

if nargin<4, tol=[]; end
if isempty(tol), tol=1e-3; end

intf=fun_sp(0,za,zb,tol);
N=real(round(intf));
if(N~=0 && ~isnan(N))
    domain=repmat(struct('za',0.0,'zb',0.0,'N',0),N+1,1);

    domain(1).za=za;
    domain(1).zb=zb;
    domain(1).N=N;

    if(N>M)
        nonzeroindex=1;
    end
    emptyindex=2:(N+1);
    while(max([domain.N])>M)
%     while(max([domain(nonzeroindex).N])>M)

        jd=nonzeroindex(1);

        tmpza=domain(jd).za;
        tmpzb=domain(jd).zb;
        tmpzc=0.49*domain(jd).za+0.51*domain(jd).zb;

        emptyindex=[jd,emptyindex];
        nonzeroindex(1)=[];
        domain(jd).N=0;

        % divide to four subdomains
        ind=emptyindex(1);
        domain(ind).za=tmpza;
        domain(ind).zb=tmpzc;
        domain(ind).N=round(fun_sp(0,domain(ind).za,domain(ind).zb,tol));
        if(domain(ind).N>0) % domain(ind).N==0 or NaN, or not interge
            emptyindex(1)=[];
        end
        if(domain(ind).N>M)
            nonzeroindex=[nonzeroindex,ind];
        end
        ind=emptyindex(1);
        domain(ind).za=real(tmpzc)+1i*imag(tmpza);
        domain(ind).zb=real(tmpzb)+1i*imag(tmpzc);
        domain(ind).N=round(fun_sp(0,domain(ind).za,domain(ind).zb,tol));
        if(domain(ind).N>0)
            emptyindex(1)=[];
        end
        if(domain(ind).N>M)
            nonzeroindex=[nonzeroindex,ind];
        end
        ind=emptyindex(1);
        domain(ind).za=real(tmpza)+1i*imag(tmpzc);
        domain(ind).zb=real(tmpzc)+1i*imag(tmpzb);
        domain(ind).N=round(fun_sp(0,domain(ind).za,domain(ind).zb,tol));
        if(domain(ind).N>0)
            emptyindex(1)=[];
        end
        if(domain(ind).N>M)
            nonzeroindex=[nonzeroindex,ind];
        end
        ind=emptyindex(1);
        domain(ind).za=tmpzc;
        domain(ind).zb=tmpzb;
        domain(ind).N=round(fun_sp(0,domain(ind).za,domain(ind).zb,tol));
        if(domain(ind).N>0)
            emptyindex(1)=[];
        end
        if(domain(ind).N>M)
            nonzeroindex=[nonzeroindex,ind];
        end
    end
    domain(emptyindex)=[];

else
    domain.za=za;
    domain.zb=za;
    domain.N=0;
end
