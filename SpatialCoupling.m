function Diff = SpatialCoupling(W,E,I,Coords)

   
    ConnMat = W;
    Ex = logical(E);
    In = logical(I);
    Dist = bsxfun(@minus,Coords,Coords');
    ProjEx = ConnMat(:,Ex); 
    ProjIn = ConnMat(:,In);
    DistEx = Dist(:,Ex);
    DistIn = Dist(:,In);
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

    Diff = abs(MExZ.*CExZ)-abs(MInZ.*CInZ);

end