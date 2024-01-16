function [SymW,BinSums] = SymmetricGains(W,Coords)

    f_exc = 0.5;
    W = W;
    L = round(range(Coords),-2);
    edges = [-L:10:L];    
    Dist = bsxfun(@minus,Coords,Coords');
    Targets = W(:,:) ~= 0;
    ProjDist = Dist.*Targets;

    BinSums = zeros(100,2);

    for ii = 1:100
        clearvars -except f_exc W L edges Dist Targets ProjDist ii Coords BinSums

        CurrNegBin = [edges(ii) edges(ii+1)];
        NegSynapses = (ProjDist >= CurrNegBin(1)) & (ProjDist < CurrNegBin(2));
        [row,col] = find(NegSynapses);
        NegLinIx= sub2ind(size(ProjDist), row, col);
        NegSynVals = W(NegLinIx);
        NegSumVals = sum(NegSynVals);

        CurrPosBin = -[edges(ii) edges(ii+1)];
        PosSynapses = (ProjDist <= CurrPosBin(1)) & (ProjDist > CurrPosBin(2));
        [row,col] = find(PosSynapses);
        PosLinIx= sub2ind(size(ProjDist), row, col);
        PosSynVals = W(PosLinIx);
        PosSumVals = sum(PosSynVals);

        if isempty(PosLinIx)
                W(NegLinIx) = 0;

        elseif isempty(NegLinIx)
                W(PosLinIx) = 0;

        elseif sign(NegSumVals) == sign(PosSumVals) && PosSumVals*NegSumVals ~= 0

            MeanSum = mean([NegSumVals PosSumVals]);
            NegScaling = MeanSum/NegSumVals;
            PosScaling = MeanSum/PosSumVals;
  
            W(NegLinIx) = W(NegLinIx)*NegScaling;
            W(PosLinIx) = W(PosLinIx)*PosScaling;


        elseif sign(NegSumVals) ~= sign(PosSumVals) && PosSumVals*NegSumVals ~= 0
            
            ExcSynIx = NegSynVals > 0;
            InhSynIx = NegSynVals < 0;

            if ~any(ExcSynIx) || ~any(InhSynIx)

                W(NegLinIx) = 0;
                W(PosLinIx) = 0;

            else

                %% Balance the (-) positioned bin
                ExcSum = sum(NegSynVals(ExcSynIx));
                InhSum = sum(NegSynVals(InhSynIx));
    
                ExcInhMean = mean(abs([InhSum ExcSum]));
    
                ExcScaler = ExcInhMean/ExcSum;
                InhScaler = ExcInhMean/InhSum;
    
                ExcLinIx = NegLinIx.*ExcSynIx;
                InhLinIx = NegLinIx.*InhSynIx;
    
                NonZero = ExcLinIx~= 0;
                ExcLinIx = ExcLinIx(NonZero);
    
                NonZero = InhLinIx~= 0;
                InhLinIx = InhLinIx(NonZero);
    
    
                W(ExcLinIx) = W(ExcLinIx)*ExcScaler;
                W(InhLinIx) = W(InhLinIx)*InhScaler; 
                W(InhLinIx) = -W(InhLinIx);
                         
    
                ExcSynIx = PosSynVals > 0;
                InhSynIx = PosSynVals < 0;

                if ~any(ExcSynIx) || ~any(InhSynIx)

                    W(NegLinIx) = 0;
                    W(PosLinIx) = 0;

                else
                %% Balance the (+) positioned bin 

                    ExcSum = sum(PosSynVals(ExcSynIx));
                    InhSum = sum(PosSynVals(InhSynIx));
        
                    ExcInhMean = mean(abs([InhSum ExcSum]));
        
                    ExcScaler = ExcInhMean/ExcSum;
                    InhScaler = ExcInhMean/InhSum;
        
                    ExcLinIx = PosLinIx.*ExcSynIx;
                    InhLinIx = PosLinIx.*InhSynIx;
        
                    NonZero = ExcLinIx~= 0;
                    ExcLinIx = ExcLinIx(NonZero);
        
                    NonZero = InhLinIx~= 0;
                    InhLinIx = InhLinIx(NonZero);
    
    
                    W(ExcLinIx) = W(ExcLinIx)*ExcScaler;
                    W(InhLinIx) = W(InhLinIx)*InhScaler; 
                    W(InhLinIx) = -W(InhLinIx);   
                end
            end
            
        end

        BinSums(ii,:) = [sum(W(NegLinIx)) sum(W(PosLinIx))];
    end

    SymW = W;
    rad = range(real(eig(SymW)));
    SymW = SymW/(rad*0.5);
%   SymW = BalanceConnectivity(SymW);
    Rates = SimulateNetwork(SymW,20000);
    Rates = Rates(10001:end,:);
    CouplingKernel(SymW,f_exc,Coords,Rates);

end

