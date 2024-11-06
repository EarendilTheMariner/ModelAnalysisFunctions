function obj = Instantiate(obj)
    gain = 1;
    Var = 0.1;

    n = obj.Length;
    Dist = -bsxfun(@minus,obj.Coordinates,obj.Coordinates');
    L = range(obj.Coordinates,2);
    edges = -L:10:L;

%% Subpopulation Indices
    IndxTemp = zeros(1,20);

    sp1 = IndxTemp;
    sp1([1, 11]) = 1;
    sp1 = logical(sp1);
    sp1 = repmat(sp1',n/20,1);

    sp5 = IndxTemp;
    sp5([2, 12]) = 1;
    sp5 = logical(sp5);
    sp5 = repmat(sp5',n/20,1);

    sp4 = IndxTemp;
    sp4([5,6,7,13,14]) = 1;
    sp4 = logical(sp4);
    sp4 = repmat(sp4',n/20,1);

    sp2 = IndxTemp;
    sp2([3,4,15,16,17]) = 1;
    sp2 = logical(sp2);
    sp2 = repmat(sp2',n/20,1);

    sp3 = IndxTemp;
    sp3([8,9,10,18,19,20]) = 1;
    sp3 = logical(sp3);
    sp3 = repmat(sp3',n/20,1);

    obj.Types(:,1) = logical(sp1);
    obj.Types(:,2) = logical(sp2);
    obj.Types(:,3) = logical(sp3);
    obj.Types(:,4) = logical(sp4);
    obj.Types(:,5) = logical(sp5);


%% Distances from subpopulations

    Dist1 = round(Dist(:,sp1),0);
    Dist2 = round(Dist(:,sp2),0);
    Dist3 = round(Dist(:,sp3),0);
    Dist4 = round(Dist(:,sp4),0);
    Dist5 = round(Dist(:,sp5),0);

    Y1 = discretize(Dist1,edges);
    Y2 = discretize(Dist2,edges);
    Y3 = discretize(Dist3,edges);
    Y4 = discretize(Dist4,edges);
    Y5 = discretize(Dist5,edges);

%%  Projection biases (Gaussians)

    DistDomain = 1:length(edges)-1;

    Sp1ProjIx = DistDomain(int32(length(DistDomain)*obj.PopParams(1)));
    Sp2ProjIx = DistDomain(int32(length(DistDomain)*obj.PopParams(2)));
    Sp3ProjIx = DistDomain(int32(length(DistDomain)*obj.PopParams(3)));
    Sp4ProjIx = DistDomain(int32(length(DistDomain)*obj.PopParams(4)));
    Sp5ProjIx = DistDomain(int32(length(DistDomain)*obj.PopParams(5)));

    sigma = length(DistDomain)/60; % Feel free to play with other Gaussian widths for the projections 
    CnSigma = sigma*(3/2);
    PDF1 = exp(-(DistDomain - Sp1ProjIx).^2 / (2 * sigma^2));
    PDF1 = [0 PDF1(1:end-1)] + PDF1;  % I do this step to improve the discrete sampling, by making the PDF symmetric around the specified mean projection distance
    PDF1 = PDF1/sum(PDF1);

    PDF2 = exp(-(DistDomain - Sp2ProjIx).^2 / (2 * sigma^2));
    PDF2 = [0 PDF2(1:end-1)] + PDF2;
    PDF2 = PDF2/sum(PDF2);

    PDF3 = exp(-(DistDomain - Sp3ProjIx).^2 / (2 * CnSigma^2));
    PDF3 = [0 PDF3(1:end-1)] + PDF3;
    PDF3 = PDF3/sum(PDF3);

    PDF4 = exp(-(DistDomain - Sp4ProjIx).^2 / (2 * sigma^2));
    PDF4 = [0 PDF4(1:end-1)] + PDF4;
    PDF4 = PDF4/sum(PDF4);

    PDF5 = exp(-(DistDomain - Sp5ProjIx).^2 / (2 * sigma^2));
    PDF5 = [0 PDF5(1:end-1)] + PDF5;
    PDF5 = PDF5/sum(PDF5);

%% Connection propabilies

    Prob1 = 15*PDF1(int32(Y1));
    Prob1(Prob1 > 1) = 1;

    Prob2 = 15*PDF2(int32(Y2));
    Prob2(Prob2 > 1) = 1;

    Prob3 = 15*PDF3(int32(Y3));
    Prob3(Prob3 > 1) = 1;

    Prob4 = 15*PDF4(int32(Y4));
    Prob4(Prob4 > 1) = 1;

    Prob5 = 15*PDF5(int32(Y5));
    Prob5(Prob5 > 1) = 1;
    

%%  SYNAPSIFY! (Sample connections from probs)

    sp1Cons = binornd(1,Prob1);

    sp2Cons = binornd(1,Prob2);

    sp3Cons = binornd(1,Prob3);

    sp4Cons = binornd(1,Prob4);

    sp5Cons = binornd(1,Prob5);

    GainMat = normrnd(gain,Var,size(sp1Cons));
    sp1Weights = sp1Cons.*abs(GainMat);


    GainMat = normrnd(gain,Var,size(sp3Cons));
    sp3Weights = sp3Cons.*abs(GainMat);


    GainMat = normrnd(gain,Var,size(sp5Cons));
    sp5Weights = sp5Cons.*abs(GainMat);


    GainMat = normrnd(gain,Var,size(sp2Cons));
    sp2Weights = -sp2Cons.*abs(GainMat);

    GainMat = normrnd(gain,Var,size(sp4Cons));
    sp4Weights = -sp4Cons.*abs(GainMat);


%% Full matrix and balancing

    W = zeros(size(Dist));

    W(:,sp1) = sp1Weights;
    W(:,sp2) = sp2Weights;
    W(:,sp3) = sp3Weights;
    W(:,sp4) = sp4Weights;
    W(:,sp5) = sp5Weights;

    [W,~] = BalanceNormalize(W);
    obj.ConnMat = W;
    obj.PDFs = [PDF1; PDF2; PDF3; PDF4; PDF5]';

    obj.E = logical(sp1 + sp3 + sp5);
    obj.I = logical(sp2 + sp4);

    return
end

