function [W,Coords,Diff,E,I] = BalancedRecStructured(Size,Length,Projectome)
    rng("shuffle")
    gain = 1;
    Var = 0.1;

    n = Size;
    W = zeros(n,n);
    Coords = linspace(0,Length,n);
    Dist = -bsxfun(@minus,Coords,Coords');
    L = range(Coords,2);
    edges = -L:10:L;
    MeanDiff = Projectome;

    E = zeros(1,Size);
    I = zeros(1,Size);
    for ii = 1:2:Size
        E(ii) = 1;
        I(ii+1) = 1;
    end
    

    DistI = round(Dist(:,logical(I)),0);
    DistE = round(Dist(:,logical(E)),0);

    YE = discretize(DistE,edges);
    YI = discretize(DistI,edges);

    % Construct PDFs for populations

    LengthDiff = (length(edges)-1)-length(MeanDiff);

    Padding = zeros(LengthDiff/2,1);
    MeanDiff = [Padding; MeanDiff; Padding];


    PosIx = MeanDiff > 0;
    NegIx = MeanDiff < 0;

    ExcPDF = zeros(length(edges)-1);
    ExcPDF(PosIx) = MeanDiff(PosIx);
    ExcPDF = ExcPDF/sum(ExcPDF);

    InhPDF = zeros(length(edges)-1);
    InhPDF(NegIx) = abs(MeanDiff(NegIx));
    InhPDF = InhPDF/sum(InhPDF);

    % Construct probability space

    ProbI = 14.5*InhPDF(YI);
    ProbE = 14.5*ExcPDF(YE);

    ProbI(ProbI>1) = 1;
    ProbE(ProbE>1) = 1;

    % Sample connections 

    InhCons = binornd(1,ProbI);
    ExcCons = binornd(1,ProbE);

    GainMatE = normrnd(gain,Var,size(ExcCons));   
    GainMatI = normrnd(gain,Var,size(InhCons)); 


    InhCons = -InhCons.*abs(GainMatI);
    ExcCons = ExcCons.*abs(GainMatE);

    W = zeros(size(Dist));
    
    W(:,logical(I)) = InhCons;
    W(:,logical(E)) = ExcCons;

    W = BalanceNormalize(W);

    Diff = SpatialCoupling(W,E,I,Coords,false);
    
   % Rates = SimulateNetwork(W,20000);
   % Rates = Rates(10001:end,:);

%    [brevr,bixev] = sort(real(Eigenvalues),'descend');
%    bievr = imag(Eigenvalues(bixev));

%    [~,eix] = max(real(Eigenvalues));
%    DominantMode = Eigenvalues(eix);

%    set(gcf, 'WindowState', 'maximized');
%    drawnow;  % Ensure the plot is fully rendered
%    save2pdf(fig,['./'],['Spectrum_2'],'-dpdf');
%    close all


end

