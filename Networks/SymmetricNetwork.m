function [W,Coords,Diff,PopMat,E,I,PDFs] = SymmetricNetwork(Size,Length,MeanDiff,AsymPopLoc)

    gain = 1;
    Var = 0.1;

    n = Size;
    W = zeros(n,n);
    Coords = linspace(0,Length,n);
    Dist = -bsxfun(@minus,Coords,Coords');
    L = range(Coords,2);
    edges = -L:10:L;

%% Subpopulation Indices

    IndxTemp = zeros(1,20);

    sp1 = IndxTemp;
    sp1([9, 19]) = 1;
    sp1 = logical(sp1);
    sp1 = repmat(sp1',n/20,1);

    sp5 = IndxTemp;
    sp5([1, 11]) = 1;
    sp5 = logical(sp5);
    sp5 = repmat(sp5',n/20,1);

    sp4 = IndxTemp;
    sp4([2,7,8,12,17]) = 1;
    sp4 = logical(sp4);
    sp4 = repmat(sp4',n/20,1);

    sp2 = IndxTemp;
    sp2([3,10,13,18,20]) = 1;
    sp2 = logical(sp2);
    sp2 = repmat(sp2',n/20,1);

    sp3 = IndxTemp;
    sp3([4,6,14,16]) = 1;
    sp3 = logical(sp3);
    sp3 = repmat(sp3',n/20,1);

    sp6 = IndxTemp;
    sp6([5,15]) = 1;
    sp6 = logical(sp6);
    sp6 = repmat(sp6',n/20,1);

    PopMat = zeros(n,6);
    PopMat(:,1) = sp1;
    PopMat(:,2) = sp2;
    PopMat(:,3) = sp3;
    PopMat(:,4) = sp4;
    PopMat(:,5) = sp5;
    PopMat(:,6) = sp6;


%% Distances from subpopulations

    Dist1 = round(Dist(:,sp1),0);
    Dist2 = round(Dist(:,sp2),0);
    Dist3 = round(Dist(:,sp3),0);
    Dist4 = round(Dist(:,sp4),0);
    Dist5 = round(Dist(:,sp5),0);
    Dist6 = round(Dist(:,sp6),0);

    Y1 = discretize(Dist1,edges);
    Y2 = discretize(Dist2,edges);
    Y3 = discretize(Dist3,edges);
    Y4 = discretize(Dist4,edges);
    Y5 = discretize(Dist5,edges);
    Y6 = discretize(Dist6,edges);

%%  Subpopulation projectome indices for PDFs

    LengthDiff = (length(edges)-1)-length(MeanDiff);

    Padding = zeros(LengthDiff/2,1);
    

    PDF1Ix = zeros(length(MeanDiff),1);
    PDF1Ix(1:15) = 1;
    PDF1Ix = [Padding; PDF1Ix; Padding];
    PDF1Ix = logical(PDF1Ix);
    

    PDF2Ix = zeros(length(MeanDiff),1);
    PDF2Ix(16:78) = 1;
    PDF2Ix = [Padding; PDF2Ix; Padding];
    PDF2Ix = logical(PDF2Ix);


    PDF3Ix = zeros(length(MeanDiff),1);
    PDF3Ix(79:122) = 1;
    PDF3Ix = [Padding; PDF3Ix; Padding];
    PDF3Ix = logical(PDF3Ix);


    PDF4Ix = zeros(length(MeanDiff),1);
    PDF4Ix(123:185) = 1;
    PDF4Ix = [Padding; PDF4Ix; Padding];
    PDF4Ix = logical(PDF4Ix);


    PDF5Ix = zeros(length(MeanDiff),1);
    PDF5Ix(186:end) = 1; 
    PDF5Ix = [Padding; PDF5Ix; Padding];
    PDF5Ix = logical(PDF5Ix);


    MeanDiff = [Padding; MeanDiff; Padding];


%% Spatial PDFs for subpopulations


    PDF1 = MeanDiff;
    PDF1(~PDF1Ix) = false;
    Domain = 1:length(MeanDiff);
    Domain = Domain';
    GaussFit = fit(Domain,PDF1,'gauss1');
    PDF1 = GaussFit(Domain);
    PDF1 = PDF1/sum(PDF1);

    PDF2 = MeanDiff;
    PDF2(~PDF2Ix) = false;
    Domain = 1:length(MeanDiff);
    Domain = Domain';
    GaussFit = fit(Domain,PDF2,'gauss1');
    PDF2 = GaussFit(Domain);
    PDF2 = PDF2/sum(PDF2);

    PDF3 = MeanDiff;
    PDF3(~PDF3Ix) = false;
    Domain = 1:length(MeanDiff);
    Domain = Domain';
    GaussFit = fit(Domain,PDF3,'gauss1');
    PDF3 = GaussFit(Domain);
    PDF3 = PDF3/sum(PDF3);

    PDF4 = MeanDiff;
    PDF4(~PDF4Ix) = false;
    Domain = 1:length(MeanDiff);
    Domain = Domain';
    GaussFit = fit(Domain,PDF4,'gauss1');
    PDF4 = GaussFit(Domain);
    PDF4 = PDF4/sum(PDF4);

    PDF5 = MeanDiff;
    PDF5(~PDF5Ix) = false;
    Domain = 1:length(MeanDiff);
    Domain = Domain';
    GaussFit = fit(Domain,PDF5,'gauss1');
    PDF5 = GaussFit(Domain);
    PDF5 = PDF5/sum(PDF5);

    Domain = 1:length(MeanDiff);
    Domain = Domain';
    y = zeros(1,length(MeanDiff));
    y(AsymPopLoc:AsymPopLoc+14) = MeanDiff(1:15);
    GaussFit = fit(Domain,y','gauss1');
    PDF6 = GaussFit(Domain);
    PDF6 = PDF6/sum(PDF6);

    PDFs = [PDF1'; PDF2'; PDF3'; PDF4'; PDF5'; PDF6'];
    PDFs = PDFs';


%% Connection propabilies

    Prob1 = 20*PDF1(Y1);
    Prob1(Prob1 > 1) = 1;
  %  Prob1(sp4,:) = Prob1(sp4,:)*0.1;

    Prob2 = 20*PDF2(Y2);
    Prob2(Prob2 > 1) = 1;
%    Prob2(sp3,:) = Prob2(sp3,:)*0.75;
%    Prob2(sp2,:) = Prob2(sp2,:)*0.1;

    Prob3 = 20*PDF3(Y3);
    Prob3(Prob3 > 1) = 1;
%    Prob3(sp4,:) = Prob3(sp4,:)*0.75;
%    Prob3(sp4,:) = Prob3(sp4,:)*0.1;

    Prob4 = 20*PDF4(Y4);
    Prob4(Prob4 > 1) = 1;
%    Prob4(sp4,:) = Prob4(sp4,:)*0.1;

    Prob5 = 20*PDF5(Y5);
    Prob5(Prob5 > 1) = 1;
 %   Prob5(sp4,:) = Prob5(sp4,:)*0.1;

    Prob6 = 20*PDF6(Y6);
    Prob6(Prob6 > 1) = 1;


    

%%  SYNAPSIFY! (Sample connections from probs)

    sp1Cons = binornd(1,Prob1);

    sp2Cons = binornd(1,Prob2);

    sp3Cons = binornd(1,Prob3);

    sp4Cons = binornd(1,Prob4);

    sp5Cons = binornd(1,Prob5);

    sp6Cons = binornd(1,Prob6);

    GainMat = normrnd(gain,Var,size(sp1Cons));
    sp1Weights = sp1Cons.*abs(GainMat);


    GainMat = normrnd(gain,Var,size(sp3Cons));
    sp3Weights = sp3Cons.*abs(GainMat);


    GainMat = normrnd(gain,Var,size(sp5Cons));
    sp5Weights = sp5Cons.*abs(GainMat);

    GainMat = normrnd(gain,Var,size(sp6Cons));
    sp6Weights = -sp6Cons.*abs(GainMat);


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
    W(:,sp6) = sp6Weights;

    [W,~] = BalanceNormalize(W);

    E = sp1 + sp3 + sp5;
    I = sp2 + sp4 + sp6;

    Diff = SpatialCoupling(W,E,I,Coords,true);


end



