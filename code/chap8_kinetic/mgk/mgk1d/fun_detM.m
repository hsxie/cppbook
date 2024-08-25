% Hua-sheng XIE, FSC-PKU, huashengxie@gmail.com, 2017-02-07 18:06
% 2017-02-12 19:47 remove J0~exp(-x^2/4) and use besselj directly
function var=fun_detM(w)
% tt0 is theta', tt is theta

global tau kt q s epsn etai tk;
global dvx dvy vxx vyy dth tt;
global nL hll Mij;

tt0=tt;
Mij=zeros(nL,nL);
for vy=vyy
    for vx=vxx        
        
        wT=-kt/epsn*(1+etai*(0.5*(vx^2+vy^2)-1.5));
        wD=-kt*(vy^2/2+vx^2);
        fm=exp(-(vx^2+vy^2)/2)/(sqrt(2*pi)^3);
		
		J0=besselj(0,kt*sqrt(1+(tt-tk).^2*s^2)*vy);

        tmp1=[];
        kappa=-2*pi*q*1i*tau*vy/vx*(w-wT)*fm;
        for jth=1:length(tt)
            phase=1i*w/vx*q*abs(tt(jth)-tt0)-1i*wD/vx*q*sign(tt(jth)-tt0).*...
                ( (1+s)*(sin(tt(jth))-sin(tt0))+s*((tk-tt(jth))*cos(tt(jth)) -...
                (tk-tt0).*cos(tt0) ));
            tmp0=J0(jth)*kappa*dth*(J0.*exp(phase))*hll.'; % integral dtheta'
            tmp1=[tmp1;tmp0];
        end
        tmp=hll*tmp1*dth; % integral dtheta

        Mij=Mij+tmp;
    end % integral dvpara
end % integral dvperp
Mij=(1+tau)*eye(nL)-Mij*dvx*dvy;
var=det(Mij);
