function obj = Instantiate2(obj)
    gain = 1;
    Var = 0.1;

    n = obj.Length;
    Dist = -bsxfun(@minus,obj.Coordinates,obj.Coordinates');
    L = range(obj.Coordinates,2);
    edges = -L:10:L;

%% Subpopulation Indices
    IndxTemp = zeros(1,10);

    sp2 = IndxTemp;
    sp2([1,4,8]) = 1;
    sp2 = logical(sp2);
    sp2 = repmat(sp2',n/10,1);

    sp3 = IndxTemp;
    sp3([3,5,7,9]) = 1;
    sp3 = logical(sp3);
    sp3 = repmat(sp3',n/10,1);

    sp4 = IndxTemp;
    sp4([2,6,10]) = 1;
    sp4 = logical(sp4);
    sp4 = repmat(sp4',n/10,1);

    obj.Types(:,2) = logical(sp2);
    obj.Types(:,3) = logical(sp3);
    obj.Types(:,4) = logical(sp4);


%% Distances from subpopulations

    Dist2 = round(Dist(:,sp2),0);
    Dist3 = round(Dist(:,sp3),0);
    Dist4 = round(Dist(:,sp4),0);

    Y2 = discretize(Dist2,edges);
    Y3 = discretize(Dist3,edges);
    Y4 = discretize(Dist4,edges);

%%  Projection biases (Gaussians)

    DistDomain = 1:length(edges)-1;

    Sp2ProjIx = DistDomain(int32(length(DistDomain)*obj.PopParams(1)));
    Sp3ProjIx = DistDomain(int32(length(DistDomain)*obj.PopParams(2)));
    Sp4ProjIx = DistDomain(int32(length(DistDomain)*obj.PopParams(3)));

    sigma = length(DistDomain)/60; % Feel free to play with other Gaussian widths for the projections 
    CnSigma = sigma*(6/2);

    PDF2 = exp(-(DistDomain - Sp2ProjIx).^2 / (2 * sigma^2));
    PDF2 = [0 PDF2(1:end-1)] + PDF2;
    PDF2 = PDF2/sum(PDF2);

    PDF3 = exp(-(DistDomain - Sp3ProjIx).^2 / (2 * CnSigma^2));
    PDF3 = [0 PDF3(1:end-1)] + PDF3;
    PDF3 = PDF3/sum(PDF3);

    PDF4 = exp(-(DistDomain - Sp4ProjIx).^2 / (2 * sigma^2));
    PDF4 = [0 PDF4(1:end-1)] + PDF4;
    PDF4 = PDF4/sum(PDF4);


%% Connection propabilies

    Prob2 = 15*PDF2(int32(Y2));
    Prob2(Prob2 > 1) = 1;

    Prob3 = 15*PDF3(int32(Y3));
    Prob3(Prob3 > 1) = 1;

    Prob4 = 15*PDF4(int32(Y4));
    Prob4(Prob4 > 1) = 1;


%%  SYNAPSIFY! (Sample connections from probs)


    sp2Cons = binornd(1,Prob2);

    sp3Cons = binornd(1,Prob3);

    sp4Cons = binornd(1,Prob4);

    GainMat = normrnd(gain,Var,size(sp2Cons));
    sp2Weights = -sp2Cons.*abs(GainMat);

    GainMat = normrnd(gain,Var,size(sp3Cons));
    sp3Weights = sp3Cons.*abs(GainMat);

    GainMat = normrnd(gain,Var,size(sp4Cons));
    sp4Weights = -sp4Cons.*abs(GainMat);


%% Full matrix and balancing

    W = zeros(size(Dist));

    W(:,sp2) = sp2Weights;
    W(:,sp3) = sp3Weights;
    W(:,sp4) = sp4Weights;

    [W,~] = BalanceNormalize(W);
    obj.ConnMat = W;
    obj.PDFs = [PDF2; PDF3; PDF4]';

    obj.E = logical(sp3);
    obj.I = logical(sp2 + sp4);

    return
end