function [Connectivity,Sparsity] = RescaleConnectivity(Connectivity,varargin)
    % Rescale Connectivity
    r = (numel(Connectivity(Connectivity>0))/numel(Connectivity(Connectivity<0)));
    Sparsity = nnz(~Connectivity)/numel(Connectivity);
    Connectivity = Connectivity./(sqrt(Sparsity*(1.-Sparsity)*(1.+(r^2))/2)*sqrt(length(Connectivity)));
end
