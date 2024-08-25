% Hua-sheng XIE, 2016-08-30 19;24
function Ep = Eptoc(Epring,np,nring)
% np*nring particles E to np guiding center E
Ep=zeros(np,1);

for ip = 1:np
    ind=nring*(ip-1);
    for ir=1:nring
        Ep(ip)=Ep(ip)+Epring(ind+ir);
    end
end
Ep=Ep/nring;