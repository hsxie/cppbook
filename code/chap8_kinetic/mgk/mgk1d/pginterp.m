% Hua-sheng XIE, 2016-04-20 15;43
% 16-08-31 16;40
function pgmat = pginterp(zp,ng,dx,xmin,np)
% interpolation between partilce positions and grids
    Ip=zeros(2*np,1); Jg=zeros(2*np,1); Wpg=zeros(2*np,1);
    for ip = 1:np
        ind=2*(ip-1);
        lindex=floor((zp(ip,1)-xmin)/dx)+1; % fisrt grid xg(1)=xmin
        dxleft=zp(ip,1)-xmin-(lindex-1)*dx;
        
        % left weight
        Ip(ind+1)=ip;
        Jg(ind+1)=lindex;
        Wpg(ind+1)=1-dxleft/dx;
        
        % right weight
        Ip(ind+2)=ip;
        Jg(ind+2)=lindex*(lindex<ng)+1;% lindex=ng, last grid # ng+1 -> 1
        Wpg(ind+2)=dxleft/dx;
    end
    pgmat=sparse(Ip,Jg,Wpg,np,ng); % matrix for interpolation
end

