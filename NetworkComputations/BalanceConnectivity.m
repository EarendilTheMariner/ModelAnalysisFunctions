function [Connectivity,Sparsity] = BalanceConnectivity(Connectivity,varargin)
    % Balance Network


    Compression = false;

    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'Compression'
                Compression = varargin{ii+1};
        end
    end


    for ii = 1:size(Connectivity,1)
        balance = sum(Connectivity(ii,:));
        jjpos = find(Connectivity(ii,:)>0);
        jjneg = find(Connectivity(ii,:)<0);

        bp = sum(Connectivity(ii,jjpos))-(balance/2);
        bn = sum(Connectivity(ii,jjneg))-(balance/2);

        facp = bp/sum(Connectivity(ii,jjpos));

        facn = bn/sum(Connectivity(ii,jjneg));


        Connectivity(ii,jjpos) = Connectivity(ii,jjpos)*facp; 
        Connectivity(ii,jjneg) = Connectivity(ii,jjneg)*facn; 

        % np = numel(jjpos);
        % nn = numel(jjneg);
        % 
        % Connectivity(ii,jjpos) = Connectivity(ii,jjpos) - (balance/(2*np)); 
        % Connectivity(ii,jjneg) = Connectivity(ii,jjneg) - (balance/(2*nn)); 
    end    

    if Compression == true
        r = (numel(Connectivity(Connectivity>0))/numel(Connectivity(Connectivity<0)));
        Sparsity = nnz(~Connectivity)/numel(Connectivity);
        Connectivity = Connectivity./(sqrt(Sparsity*(1.-Sparsity)*(1.+(r^2))/2)*sqrt(length(Connectivity)));
    else
        disp("Raw Balance")
    end

end
