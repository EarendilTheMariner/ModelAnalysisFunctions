function [W,Eigenvalues] = BalanceNormalize(W)


    W = BalanceConnectivity(W);

    Eigenvalues = eig(W);

    [~,eix] = max(real(Eigenvalues));
    
    DominantMode = Eigenvalues(eix);
    
    W = W/(real(DominantMode)*0.95);
    
    Eigenvalues = eig(W);

 %   ScaleFactor = 1/(real(DominantMode)*0.95);
 %   save("GlobalScaling.mat","ScaleFactor");
    

end

