function obj = GainMod(obj,Pop,Gain,normalizeFlag)

    % Method for scaling the synaptic gain of subpopulation of interest


   

    if Pop == 2

        Pop = obj.Types(:,Pop);
        Pop = logical(Pop);
        W = obj.ConnMat;
        S1 =  sum(W(:,Pop),"all");
        S2 = sum(W(:,logical(obj.Types(:,4))),"all");
        STot = S1 + S2;
        display(S1);
        display(S2);
        display(STot);

        W(:,Pop) = W(:,Pop).*Gain;
        CounterGain = (STot-Gain*S1)/S2;    
        W(:,logical(obj.Types(:,4))) = W(:,logical(obj.Types(:,4))).*CounterGain;

     
        S1 =  sum(W(:,Pop),"all");
        S2 = sum(W(:,logical(obj.Types(:,4))),"all");
        STot = S1 + S2;
        display(S1);
        display(S2);
        display(STot);


    elseif Pop == 4

        Pop = obj.Types(:,Pop);
        Pop = logical(Pop);
        W = obj.ConnMat;
        S2 =  sum(W(:,Pop),"all");
        S1 = sum(W(:,logical(obj.Types(:,2))),"all");
        STot = S1 + S2;

        W(:,Pop) = W(:,Pop).*Gain;
        CounterGain = (STot-Gain*S2)/S1;
        W(:,logical(obj.Types(:,2))) = W(:,logical(obj.Types(:,2))).*CounterGain;

    else

        Pop = obj.Types(:,Pop);
        Pop = logical(Pop);
        W = obj.ConnMat;
        W(:,Pop) = W(:,Pop).*Gain;
    end


    if normalizeFlag

        W = BalanceNormalize(W); % Rebalances the weights after scaling and normalized to fixed radius
    else
        W = BalanceNormalize(W,'Normalize',false); % Balance with no rescaling/normalization
    end

    obj.ConnMat = W; % Re-Initializes the connectivity matrix property after modulation

    obj.EigenMode = EigenStructure(obj.ConnMat,1); % Recomputes dominant mode

    obj.Projectome = SpatialCoupling(obj.ConnMat,obj.E,obj.I,obj.Coordinates,false);

%    obj.Skewness = SkewnessFactor(obj.Projectome);  


end

