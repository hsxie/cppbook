% By Li Yue-Yan
deltatheta=pi/101*9;
theta = [-10*pi:deltatheta:10*pi];
lengththeta = length(theta);
sintheta = sin(theta);
costheta = cos(theta);
alphacontainer = [0.005:0.005:4];
lengthalphacontainer = length(alphacontainer);
resalphacontainerlength = 0;


for alpha = alphacontainer

    fprintf('\b\b\b\b\b\b\b\b\b%g%%',alpha/4*100);
    mat01 = (1+alpha^2*sintheta.^2)./deltatheta^2;
    mat02 = alpha^2*sintheta.*costheta./deltatheta;
    mat03 = alpha*(costheta-alpha*sintheta.^2);
    mat11 = -2*alpha*theta.*sintheta./deltatheta^2;
    mat12 = -alpha*(theta.*costheta+sintheta)./deltatheta;
    mat13 = alpha*theta.*sintheta;
    mat21 = theta.^2./deltatheta^2;
    mat22 = theta./deltatheta;
    
    matL = zeros(lengththeta*2-4,lengththeta*2-4);
    matR = matL;
    
    mat0 = diag(mat03-2*mat01,0)+diag(mat01(1:end-1)+mat02(1:end-1),1)+diag(mat01(2:end)-mat02(2:end),-1);
    mat1 = diag(mat13-2*mat11,0)+diag(mat11(1:end-1)+mat12(1:end-1),1)+diag(mat11(2:end)-mat12(2:end),-1);
    mat2 = diag(-2*mat21,0)     +diag(mat21(1:end-1)+mat22(1:end-1),1)+diag(mat21(2:end)-mat22(2:end),-1);
    
    mat0=mat0(2:end-1,2:end-1);
    mat1=mat1(2:end-1,2:end-1);
    mat2=mat2(2:end-1,2:end-1);
    
    %assemble matrix
    
    matL(1:lengththeta-2                , 1:lengththeta-2              ) = eye(lengththeta-2);
    matR(1:lengththeta-2                , lengththeta+1-2:lengththeta*2-4) = eye(lengththeta-2);
    matL(lengththeta+1-2:lengththeta*2-4, lengththeta+1-2:lengththeta*2-4) = mat2;
    matR(lengththeta+1-2:lengththeta*2-4, 1:lengththeta-2              ) = -mat0;
    matR(lengththeta+1-2:lengththeta*2-4, lengththeta+1-2:lengththeta*2-4) = -mat1;
    
    matR = matL\matR;
    
    eigenscontainer = eig(matR);
    realscontainer = eigenscontainer(imag(eigenscontainer)==0);
    reslength = length(realscontainer);
    resalphacontainer(resalphacontainerlength+1:resalphacontainerlength+reslength)=ones(1,reslength)*alpha;
    resscontainer(resalphacontainerlength+1:resalphacontainerlength+reslength) = realscontainer;
    resalphacontainerlength = resalphacontainerlength+reslength;
end

plotalphacontainer = resalphacontainer(resscontainer>=0 & resscontainer<=2);
plotscontainer = resscontainer(resscontainer>=0 & resscontainer<=2);
%%
h=figure('unit','normalized','Position',[0.01 0.1 0.3 0.4],...
    'DefaultAxesFontSize',15);
plot(plotalphacontainer,plotscontainer,'k.');
ylabel('s'); xlabel('\alpha');
print(gcf,'-dpng','ballooningeignmethod2.png');