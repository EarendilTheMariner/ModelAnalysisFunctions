function SF = SkewnessFactor(AsymNet,SymNet)


    % Symmetric reference projectome
    ConnMat = SymNet.ConnMat;
    Sp1Sp5 = SymNet.Types(:,1)+SymNet.Types(:,5);
    Sp1Sp5 = logical(Sp1Sp5);
    ConnMat(:,Sp1Sp5) = ConnMat(:,Sp1Sp5).*0;

    Ex = logical(SymNet.E);
    In = logical(SymNet.I);
    Dist = -bsxfun(@minus,SymNet.Coordinates,SymNet.Coordinates');
    ProjEx = ConnMat(:,Ex); 
    ProjIn = ConnMat(:,In);
    DistEx = Dist(:,Ex);
    DistIn = Dist(:,In);
    L = round(range(SymNet.Coordinates),-2);
    edges = [-L:1:L];
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

    SymProjectome = abs(MExZ.*CExZ)-abs(MInZ.*CInZ);
    SymProjectome = smoothdata(SymProjectome/norm(SymProjectome));
    Derivative = diff(SymProjectome);
    
    
    % Asymmetric projectome
    ConnMat = AsymNet.ConnMat;
    Sp1Sp5 = AsymNet.Types(:,1)+AsymNet.Types(:,5);
    Sp1Sp5 = logical(Sp1Sp5);
    ConnMat(:,Sp1Sp5) = ConnMat(:,Sp1Sp5).*0;

    Ex = logical(AsymNet.E);
    In = logical(AsymNet.I);
    Dist = -bsxfun(@minus,AsymNet.Coordinates,AsymNet.Coordinates');
    ProjEx = ConnMat(:,Ex); 
    ProjIn = ConnMat(:,In);
    DistEx = Dist(:,Ex);
    DistIn = Dist(:,In);
    L = round(range(AsymNet.Coordinates),-2);
    edges = [-L:1:L];
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

    AsymProjectome = abs(MExZ.*CExZ)-abs(MInZ.*CInZ);
    AsymProjectome = smoothdata(AsymProjectome/norm(AsymProjectome));
    SF = -dot(Derivative,AsymProjectome(1:end-1));
   

end

 