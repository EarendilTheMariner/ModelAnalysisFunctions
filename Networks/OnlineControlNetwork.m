function [W,Coords,Diff,PopMat,E,I] = OnlineControlNetwork(Size,Length,MeanDiff)

    gain = 1;
    Var = 0.1;

    n = Size;
    W = zeros(n,n);
    Coords = linspace(0,Length,n);
    Dist = bsxfun(@minus,Coords,Coords');
    L = range(Coords,2);
    edges = -L:10:L;

%% Subpopulation Indices

    indx = 1:5:1000;

    sp1 = zeros(1,n);
    sp1(indx) = true;
    sp1 = logical(sp1);

    sp2 = zeros(1,n);
    sp2(indx+1) = true;
    sp2 = logical(sp2);

    sp3 = zeros(1,n);
    sp3(indx+2) = true;
    sp3 = logical(sp3);

    sp4 = zeros(1,n);
    sp4(indx+3) = true;
    sp4 = logical(sp4);

    sp5 = zeros(1,n);
    sp5(indx+4) = true;
    sp5 = logical(sp5);

    PopMat = zeros(L,5);
    PopMat(:,1) = sp1;
    PopMat(:,2) = sp2;
    PopMat(:,3) = sp3;
    PopMat(:,4) = sp4;
    PopMat(:,5) = sp5;

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

%%  Subpopulation projectome indices for PDFs

    PDF1Ix = MeanDiff > 0;
    PDF1Ix(round(length(MeanDiff)/4):end) = 0;

    PDF2Ix = MeanDiff < 0;
    PDF2Ix(round(length(MeanDiff)/2):end) = 0;

    PDF3Ix = MeanDiff > 0;
    PDF3Ix(1:round(length(MeanDiff)/4)) = 0;
    PDF3Ix(3*round(length(MeanDiff)/4):end) = 0;

    PDF4Ix = MeanDiff < 0;
    PDF4Ix(1:round(length(MeanDiff)/2)) = 0;

    PDF5Ix = MeanDiff > 0;
    PDF5Ix(1:3*round(length(MeanDiff)/4)) = 0;


%% Spatial PDFs for subpopulations


    PDF1 = zeros(1,length(edges)-1);
    PDF1(PDF1Ix) = MeanDiff(PDF1Ix);
    PDF1 = PDF1/sum(PDF1);


    PDF2 = zeros(1,length(edges)-1);
    PDF2(PDF2Ix) = MeanDiff(PDF2Ix);
    PDF2 = PDF2/sum(PDF2);

    PDF3 = zeros(1,length(edges)-1);
    PDF3(PDF3Ix) = MeanDiff(PDF3Ix);
    PDF3 = PDF3/sum(PDF3);
    

    PDF4 = zeros(1,length(edges)-1);
    PDF4(PDF4Ix) = MeanDiff(PDF4Ix);
    PDF4 = PDF4/sum(PDF4);

    PDF5 = zeros(1,length(edges)-1);
    PDF5(PDF5Ix) = MeanDiff(PDF5Ix);
    PDF5 = PDF5/sum(PDF5);



%% Connection propabilies

    Prob1 = 25*PDF1(Y1);
    Prob1(Prob1 > 1) = 1;

    Prob2 = 25*PDF2(Y2);
    Prob2(Prob2 > 1) = 1;
    Prob2(sp3,:) = Prob2(sp3,:)*0.75;

    Prob3 = 25*PDF3(Y3);
    Prob3(Prob3 > 1) = 1;

    Prob4 = 25*PDF4(Y4);
    Prob4(Prob4 > 1) = 1;
    Prob4(sp3,:) = Prob4(sp3,:)*0.75;

    Prob5 = 25*PDF5(Y5);
    Prob5(Prob5 > 1) = 1;
    

%%  SYNAPSIFY! (Sample connections from probs)

    sp1Cons = binornd(1,Prob1);

    sp2Cons = binornd(1,Prob2);

    sp3Cons = binornd(1,Prob3);

    sp4Cons = binornd(1,Prob4);

    sp5Cons = binornd(1,Prob5);

    GainMat = normrnd(gain,Var,size(sp1Cons));

    sp1Weights = sp1Cons.*abs(GainMat);
    sp3Weights = sp3Cons.*abs(GainMat);
    sp5Weights = sp5Cons.*abs(GainMat);

    sp2Weights = -sp2Cons.*abs(GainMat);
    sp4Weights = -sp4Cons.*abs(GainMat);


%% Full matrix and balancing

    W = zeros(size(Dist));

    W(:,sp1) = sp1Weights;
    W(:,sp2) = sp2Weights;
    W(:,sp3) = sp3Weights;
    W(:,sp4) = sp4Weights;
    W(:,sp5) = sp5Weights;

    [W,~] = BalanceNormalize(W);

    E = sp1 + sp3 + sp5;
    I = sp2 + sp4;

    Diff = SpatialCoupling(W,E,I,Coords);


end















    


    







