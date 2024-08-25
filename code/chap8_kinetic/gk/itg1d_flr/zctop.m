% Hua-sheng XIE, 2016-08-30 16;34
% 16-10-02 14:41
function zpring = zctop(zp,np,nring,L)
% np guiding center to np*nring particles
zpring=zeros(np*nring,2);

% wgt=[1/3,1/6,1/6,1/3]; % not [1/4,1/4,1/4,1/4] ref Cummings thesis p65

wgtequal=2;
if(wgtequal==1) % 2D version, not accurate for 1D
    ang=(1:1:nring).*2*pi/nring;
    wgt=0.*ang+1/nring;
	xxt=sin(ang);
elseif(wgtequal==2)
    ang=(2*(1:1:nring)-1)/nring*pi;
    ang2=(2*(1:1:nring)+1)/nring*pi;
    xxt=(sin(ang2)-sin(ang))./(ang2-ang);
    wgt=0.*ang+1/nring;
else % Cummings' nring=4
    wgt=[1/3,1/6,1/6,1/3];
    ang=[-pi/3,-pi/12,pi/12,pi/3];
	xxt=sin(ang);
end

for ip = 1:np
    ind=nring*(ip-1);
    for ir=1:nring % 16-10-02 rhoi=sqrt(2)*v_perp, not 1*v_verp!
        zpring(ind+ir,1)=zp(ip,1)+sqrt(2)*zp(ip,3)*xxt(ir); % rhoi=sqrt(2)*zp(:,3);
        zpring(ind+ir,2)=zp(ip,5)*wgt(ir); % weight
    end
end
zpring(:,1)=mod(zpring(:,1)+10*L,L);

