function [W,E] = BalanceNormalize(W,varargin)
    % Imposes row wise (detailed) balance on W and scales the weights so
    % if flag is true, scales W such that leading mode is slighlty suprastable, 

    W = BalanceConnectivity(W);
    Normalize = true;

    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'Normalize'
                Normalize = varargin{ii+1};
        end
    end
    
    if Normalize
        E = eigs(W,1);        
        W = W/(real(E)*0.95);
    end
end

