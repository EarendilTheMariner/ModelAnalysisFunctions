function [SymW,SymCoords] = SymmetricSpatialNetwork(W,Coords)
    %% 
    W = W;
    n = length(W(:,1));
    L = round(range(Coords),-2);
    edges = [-L:10:L];
    Dist = bsxfun(@minus,Coords,Coords');
    Targets = W(:,:) ~= 0;
    ProjDist = Dist.*Targets;

    for ii = 1:n
        ixPos = find(ProjDist(:,1) > 0);
        ixNeg = find(ProjDist(:,1) < 0);

        PruneIx = abs(ProjDist(:,1)) > Coords(1);

        W(PruneIx,1) = 0;



           


        if length()

            W(ixNeg,1) = 0;

            for i = 1:length(ixPos)
    
                CurrDist = ixPos(1)
                Bounds = -[ceil(CurrDist/10)*10 - 10, ceil(CurrDist/10)*10]




            




        end
    end

        


    



end 