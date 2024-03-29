function [SymW,ProjDist] = SymmetricSpatialNetwork(W,Coords)
    %% 
    f_exc = 0.5;
    gain = 1;
    var = 0.1;
    W = W;
    n = length(W(:,1));
    L = round(range(Coords),-2);
    edges = [-L:10:L];

    Dist = bsxfun(@minus,Coords,Coords');
    Targets = W(:,:) ~= 0;
    ProjDist = Dist.*Targets;


    for ii = 1:n
        
        MinDist = min(Coords(ii), (max(Coords) - Coords(ii)));
        PruneIx = ProjDist(:,ii) > MinDist | ProjDist(:,ii) < -MinDist;
        ProjDist(PruneIx,ii) = 0;
        ixPos = find(ProjDist(:,ii) > 0);
        ixNeg = find(ProjDist(:,ii) < 0);

        if max(abs(ProjDist(ixPos,ii))) > max(abs(ProjDist(ixNeg,ii)))
            ProjDist(ixNeg,ii) = 0;
        else
            ProjDist(ixPos,ii) = 0;
        end
        
        NewIx = find(ProjDist(:,ii) ~= 0);
        

        for i = 1:length(NewIx)
            CurrDist = ProjDist(NewIx(i),ii);
            Bounds = -[ceil(CurrDist/10)*10 - 10, ceil(CurrDist/10)*10];
            MaxBound = max(Bounds);
            MinBound = min(Bounds);
            SymIx = find(Dist(:,ii) <= MaxBound & Dist(:,ii) >= MinBound);

            if isempty(SymIx)
                ProjDist(NewIx(i),ii) = 0;            

            elseif length(SymIx) == 1 & ~ismember(Dist(SymIx,ii),ProjDist(:,ii))

                ProjDist(SymIx,ii) = Dist(SymIx,ii);

            elseif length(SymIx) == 1 & ismember(Dist(SymIx,ii),ProjDist(:,ii))

                ProjDist(SymIx,ii) = ProjDist(SymIx,ii);

            else

                NotInIx = find(~ismember(Dist(SymIx,ii),ProjDist(:,ii)));            
                SymIx = SymIx(NotInIx);

                if isempty(SymIx)
                    ProjDist(NewIx(i),ii) = 0;

                elseif length(SymIx) == 1
                    ProjDist(SymIx,ii) = Dist(RanSymIx,ii);
                else
                    RanSymIx = randsample(SymIx,1);
                    ProjDist(RanSymIx,ii) = Dist(RanSymIx,ii); 

                end
            end
        end
    end

    fE = f_exc*n;
    fI = (1-f_exc)*n;
    E = false(1,n);
    I = E;
    E(1:fE) = true;
    I(fE+1:end) = true;

    Connectivity = normrnd(gain,var,size(ProjDist)); 
    Connectivity(:,E) = abs(Connectivity(:,E));
    Connectivity(:,I) = -(Connectivity(:,I))*1.1;
    Filt = ProjDist ~= 0;

    SymW = Connectivity.*Filt;
    SymW = BalanceConnectivity(SymW);

end

