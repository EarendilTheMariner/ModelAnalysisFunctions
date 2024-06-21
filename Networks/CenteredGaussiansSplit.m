function [W,Coords,PopMat,E,I] = CenteredGaussiansSplit(Size,Length,Shift)
    gain = 1;
    Var = 0.1;
  
    n = Size;
    W = zeros(n,n);
    Coords = linspace(0,Length,n);
    Dist = bsxfun(@minus,Coords,Coords');
    L = range(Coords,2);
    edges = -L:10:L;

%% Subpopulation Indices

    IndxTemp = zeros(1,20);

    Exc1 = IndxTemp;
    Exc1([1, 5, 9, 13, 17]) = 1;
    Exc1 = logical(Exc1);
    Exc1 = repmat(Exc1',n/20,1);

    Exc2 = IndxTemp;
    Exc2([2, 6, 10, 14, 18]) = 1;
    Exc2 = logical(Exc2);
    Exc2 = repmat(Exc2',n/20,1);

    Inh1 = IndxTemp;
    Inh1([3, 4, 7, 8, 11, 12, 15, 16, 19, 20]) = 1;
    Inh1 = logical(Inh1);
    Inh1 = repmat(Inh1',n/20,1);


    PopMat = zeros(n,3);
    PopMat(:,1) = Exc1;
    PopMat(:,2) = Exc2;
    PopMat(:,3) = Inh1;


%% Distances from subpopulations

    Dist1 = round(Dist(:,Exc1),0);
    Dist2 = round(Dist(:,Exc2),0);
    Dist3 = round(Dist(:,Inh1),0);


    Y1 = discretize(Dist1,edges);
    Y2 = discretize(Dist2,edges);
    Y3 = discretize(Dist3,edges);


%%  Subpopulation projectome indices for PDFs

 %   LengthDiff = (length(edges)-1)-length(MeanDiff);

 %   Padding = zeros(LengthDiff/2,1);
    


%% Spatial PDFs for subpopulations


    PDF1 = normpdf(edges,Shift,50);
    PDF1 = PDF1/sum(PDF1);

    PDF2 = normpdf(edges,0,1000);
    PDF2 = PDF2/sum(PDF2);

    PDF3 = normpdf(edges,0,300);
    PDF3 = PDF3/sum(PDF3);


%% Connection propabilies

    Prob1 = 12*PDF1(Y1);
    Prob1(Prob1 > 1) = 1;
  %  Prob1(sp4,:) = Prob1(sp4,:)*0.1;

    Prob2 = 12*PDF2(Y2);
    Prob2(Prob2 > 1) = 1;
%    Prob2(sp3,:) = Prob2(sp3,:)*0.75;
%    Prob2(sp2,:) = Prob2(sp2,:)*0.1;

    Prob3 = 12*PDF3(Y3);
    Prob3(Prob3 > 1) = 1;
%    Prob3(sp4,:) = Prob3(sp4,:)*0.75;
%    Prob3(sp4,:) = Prob3(sp4,:)*0.1;


%%  SYNAPSIFY! (Sample connections from probs)

    Exc1Cons = binornd(1,Prob1);

    Exc2Cons = binornd(1,Prob2);

    Inh1Cons = binornd(1,Prob3);


    GainMat = normrnd(gain,Var,size(Exc1));
    Exc1Weights = Exc1Cons.*abs(GainMat);

    GainMat = normrnd(gain,Var,size(Exc2));
    Exc2Weights = Exc2Cons.*abs(GainMat);

    GainMat = normrnd(gain,Var,size(Inh1));
    Inh1Weights = Inh1Cons.*abs(GainMat)*(-1);



%% Full matrix and balancing

    W = zeros(size(Dist));

    W(:,Exc1) = Exc1Weights;
    W(:,Exc2) = Exc2Weights;
    W(:,Inh1) = Inh1Weights;
    [W,~] = BalanceNormalize(W);

    E = Exc1 + Exc2;
    I = Inh1;

    Diff = SpatialCoupling(W,E,I,Coords);
end

