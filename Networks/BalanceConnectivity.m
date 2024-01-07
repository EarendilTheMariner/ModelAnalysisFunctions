function [Connectivity,Sparsity] = BalanceConnectivity(Connectivity)
    % Balance Network
    for ii = 1:size(Connectivity,1)
        balance = sum(Connectivity(ii,:));
        jjpos = find(Connectivity(ii,:)>0);
        jjneg = find(Connectivity(ii,:)<0);
        np = numel(jjpos);
        nn = numel(jjneg);
        Connectivity(ii,jjpos) = Connectivity(ii,jjpos) - (balance/(2*np)); 
        Connectivity(ii,jjneg) = Connectivity(ii,jjneg) - (balance/(2*nn)); 
        %Connectivity(ii,jj) = Connectivity(ii,jj) - mean(Connectivity(ii,jj));
    end    
    r = (numel(Connectivity(Connectivity>0))/numel(Connectivity(Connectivity<0)));
    Sparsity = nnz(~Connectivity)/numel(Connectivity);
    Connectivity = Connectivity./(sqrt(Sparsity*(1.-Sparsity)*(1.+(r^2))/2)*sqrt(length(Connectivity)));
end
