function CouplingKernel(W,f_exc,Coords,Rates)

    n = length(W(:,1));
    fE = f_exc*n;
    fI = (1-f_exc)*n;
    E = false(1,n);
    I = E;
    E(1:fE) = true;
    I(fE+1:end) = true;

    screensize = get(groot,'ScreenSize');

    ConnMat = W;
    Ex = E;
    In = I;
    Dist = bsxfun(@minus,Coords,Coords');
    ProjEx = ConnMat(:,Ex); 
    ProjIn = ConnMat(:,In);
    DistEx = Dist(:,Ex);
    DistIn = Dist(:,In);
    [vec,bevr] = eig(ConnMat,'vector');
    vec = real(vec);
    [brevr,bixev] = sort(real(bevr),'descend');
    bievr = imag(bevr(bixev));
    [~,ind] = unique(brevr,'stable');
    L = round(range(Coords),-2);
    edges = [-L:10:L];
    yin = discretize(DistIn(ProjIn~= 0), edges);
    yex = discretize(DistEx(ProjEx~= 0), edges);
    [GnIn,Min,Cin] = grpstats(ProjIn(ProjIn~=0), yin,["gname","mean","numel"]);
    [GnEx,Mex,Cex] = grpstats(ProjEx(ProjEx~=0), yex,["gname","mean","numel"]);
    MInZ = zeros(length(edges)-1,1);
    MExZ = zeros(length(edges)-1,1);
    CInZ = zeros(length(edges)-1,1);
    CExZ = zeros(length(edges)-1,1);
    MInZ(str2double(GnIn)) = Min;
    MExZ(str2double(GnEx)) = Mex;
    CInZ(str2double(GnIn)) = Cin;
    CExZ(str2double(GnEx)) = Cex;
    fig = figure(Position=screensize);
    subplot (2,3,1)
    bar(edges(1:end-1), abs(MInZ.*CInZ),'FaceColor',[Colors().BergOrange],'FaceAlpha',0.5);
    hold on
    bar(edges(1:end-1), abs(MExZ.*CExZ),'FaceColor',[Colors().BergElectricBlue],'FaceAlpha',0.5);
    legend({'Inhib','Excit'});
    set(gca, 'xtick', round(edges(1:end-1),0));
    set(gca, 'xticklabels',round(edges(1:end-1),0));
    xlabel("Distance from Soma [um]");
    difference = abs(MExZ.*CExZ)-abs(MInZ.*CInZ);
    difplus = difference;
    difneg = difference;
    difplus(difplus < 0) = 0;
    difneg(difneg > 0) = 0;
    difneg = abs(difneg);
    col = [normalize(difneg,'range').*Colors().BergOrange+normalize(difplus,'range').*Colors().BergElectricBlue];
    subplot (2,3,2)
    b = bar(edges(1:end-1),abs(MExZ.*CExZ)-abs(MInZ.*CInZ),'FaceColor','flat');
    for k = 1:size(MExZ,1)
        b.CData(k,:) = col(k,:);
    end
    set(gca, 'xtick', round(edges(1:100:end-1),0));
    set(gca, 'xticklabels',round(edges(1:100:end-1),0));
    legend({'Diff Excit Inihib'});
    xlabel("Distance from Soma [um]");
    subplot (2,3,3)
    scatter(brevr,bievr,'filled','MarkerFaceColor',Colors().BergBlack,'MarkerFaceAlpha',0.5)
    vline(1)
    axis equal
    axis([-2 2 -2 2])
    subplot(2,3,4:6)
    plot(Rates(1:end,:));

end
    