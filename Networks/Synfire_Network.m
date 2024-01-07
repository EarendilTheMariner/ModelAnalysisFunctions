%%
function [Connectivity,Position,Sparsity] = Synfire_Network(varargin)
    n = 1000;
    f_exc = 0.5;
    verbose = 1;
    gainex = 1;
    gainin = 1.1;
    var = 0.1;
    length_scale_ex = 200;
    length_scale_in = 200;
    RC = 0;
    Position = nan;
    Trans = nan;
    
    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'Size'
                n = varargin{ii+1};
            case 'f_exc'
                f_exc = varargin{ii+1};
            case 'Verbose'
                verbose = varargin{ii+1};
            case 'GainEx'
                gainex = varargin{ii+1};
            case 'GainIn'
                gainin = varargin{ii+1};
            case 'LengthScaleEx'
                length_scale_ex = varargin{ii+1};
            case 'LengthScaleIn'
                length_scale_in = varargin{ii+1};
            case 'RandomComp'
                RC = varargin{ii+1};
            case 'Position'
                Position = varargin{ii+1};
            case 'Transmitter'
                Trans = varargin{ii+1};
        end 
    end
         
    if(~any(Position) & ~any(Trans))
        r = 10*rand([1 n]);
        phi = 2*pi*rand([1 n]);
        Pos = 1:2*n;
        PosX = randsample(Pos,n,1) + normrnd(0,0.2,[1 n]);
        PosY = r.*sin(phi);
        PosZ = r.*cos(phi);
        Position = [PosX' PosY' PosZ'];

        fE = f_exc*n;
        fI = (1-f_exc)*n;
       
        E = false(1,n);
        I = E;
        
        E(1:fE) = true;
        I(fE+1:end) = true;
    elseif(any(Position) & any(Trans))
        PosX = Position(:,2)';
        PosY = Position(:,1)';
        PosZ = Position(:,3)';
        
        n = length(PosX);
        f_exc = nnz(Trans>0)./numel(Trans);
          
        E = false(1,n);
        I = E;
        
        E(Trans > 0) = true;
        I(Trans < 0) = true;
    else
        error('Missing Position or Transmitter specification input')
    end

    Dist = -bsxfun(@minus,PosX,PosX');    
    DistN = zeros(size(Dist));
    DistN(:,I)= -abs(Dist(:,I)./length_scale_in);
    DistN(Dist>0 & E) = -Dist(Dist>0 & E)./length_scale_ex;
    DistN(Dist<0 & E) = Dist(Dist<0 & E)./(length_scale_ex/10);
    Prob = exp(DistN);
    
    if RC
        RandComp = zeros(size(Prob));
        ixs = randsample(1:1:numel(RandComp),round(numel(RandComp)/2),0);
        RandComp(ixs) = normrnd(0,0.5,[1 numel(ixs)]);
        Prob = Prob + RandComp;
        Prob(Prob>1) = 1;
        Prob(Prob<0) = 0;
    end
    
    %%
    Sparsifier = binornd(1,Prob); % Sparsify based on binomial proba with E = connectivity probability
    Connectivity = normrnd(gainex,var,size(Prob)); % Attribute random strength with mean gain p
    Connectivity(:,E) = abs(Connectivity(:,E));
    Connectivity(:,I) = -(Connectivity(:,I))*gainin;
    
    Connectivity = Connectivity.*Sparsifier;
    Connectivity = Connectivity - diag(diag(Connectivity)); % Remove self synapse
    Connectivity = BalanceConnectivity(Connectivity);
    Sparsity = nnz(~Connectivity)/numel(Connectivity);
    
    if verbose
        bevr = eig(Connectivity);
        [brevr,bixev] = sort(real(bevr),'descend');
        bievr = imag(bevr(bixev));
    
        fig = figure;
        scatter(brevr,bievr)
        vline(1)
        axis equal
    end
end


