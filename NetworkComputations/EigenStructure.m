function EigenStructures = EigenStructure(W,k)
%EIGENSTRUCTURE Summary of this function goes here
    [V,D] = eig(W);
    EVs = diag(D);
    [~,ind] = sort(real(EVs),'descend');
    Vs = V(:,ind);
    ModeIx = 1:2:k*2;
    EigenStructures = real(Vs(:,ModeIx));


end

