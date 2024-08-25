% Hua-sheng XIE, 2016-08-30 17;44
function zcpmat = zcpinterp(zpring,ng,dx,np,nring,L)
    np2=np*nring;
    zpring(:,1)=mod(zpring(:,1)+10*L,L);
% interpolation between partilce positions and grids
    Ip=zeros(2*np2,1); Jg=zeros(2*np2,1); Wpg=zeros(2*np2,1);
    for ip = 1:np2
        ind=2*(ip-1);
        lindex=floor(zpring(ip,1)/dx)+1; % fisrt grid xg(1)=0
%         dxleft=zp(ip,1)-lindex*dx; % wrong
        dxleft=zpring(ip,1)-(lindex-1)*dx;
        
        % left weight
        Ip(ind+1)=ip;
        Jg(ind+1)=lindex;
        Wpg(ind+1)=1-dxleft/dx;
        
        % right weight
        Ip(ind+2)=ip;
        Jg(ind+2)=lindex*(lindex<ng)+1;% lindex=ng, last grid # ng+1 -> 1
        Wpg(ind+2)=dxleft/dx;
    end
    zcpmat=sparse(Ip,Jg,Wpg,np2,ng); % matrix for interpolation
end

