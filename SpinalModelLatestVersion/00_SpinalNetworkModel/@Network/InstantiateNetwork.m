function obj = InstantiateNetwork(obj,Geometry,Parameters,varargin)
    verbose = 1;
    var = 0.1;
    RC = 0;
    SP = 0;
    Bal = 1;
    resc = 1;
    
    if(~isempty(varargin))
        for ii = 1:2:length(varargin)
            switch varargin{ii}
                case 'Verbose'
                    verbose = varargin{ii+1};
                case 'RandomComp'
                    RC = 1;
                    RCvar= varargin{ii+1};
                case 'Sparsify'
                    SP = 1;
                    sparsity = varargin{ii+1};
                case 'Balance'
                    Bal = varargin{ii+1};
                case 'Rescale'
                    resc = varargin{ii+1};
            end 
        end
    end

    PopModel = false(size(Geometry,1),1);
    for ii = 1:length(Parameters.Types)
        PopModel = PopModel | (Geometry.Type == Parameters.Types(ii));
    end 

    SC = Geometry(PopModel,:);
    BiasPop = ComputePopulationProjectionBias(SC,Parameters);
    DelaysPop = ComputePopulationDelays(SC,Parameters);
%% Generate Bias Map 
    PosX = SC.Position(:,2)';
    PosY = SC.Position(:,1)';
    PosZ = SC.Position(:,3)';

    E = (SC.Transmit > 0); 
    I = (SC.Transmit < 0); 

    CellIndex = cellfun(@(x,y) x.*(SC.Type == y),{1:length(Parameters.Types)},{Parameters.Types},'UniformOutput',false);
    CellIndex = sum(CellIndex{1,1},2);

    DistRC = -bsxfun(@minus,PosX,PosX');   
    DistML = -bsxfun(@minus,PosY,PosY');
    DistML = SC.Latera'.*DistML;
    DistDV = -bsxfun(@minus,PosZ,PosZ');
    DistCI = SC.Latera*SC.Latera';

    BiasRC = nan(size(DistRC));
    BiasML = nan(size(DistML));
    BiasDV = nan(size(DistDV));
    BiasCI = ones(size(DistCI));
    BiasLay = ones(size(DistCI));
    BiasSeg = ones(size(DistCI));
    BiasMn = ones(size(DistCI));
    BiasModule = ones(size(DistCI));


    bRC = unique(Parameters.BiasRC);
    bML = unique(Parameters.BiasML);
    bDV = unique(Parameters.BiasDV);
    bCI = unique(Parameters.BiasContraIpsi);
    bLay = Parameters.BiasLayer;
    bSeg = Parameters.BiasSegment;
    bMn = Parameters.BiasMN;

    for brc = bRC
        ioi = ismember(SC.Type,Parameters.Types(Parameters.BiasRC == brc))';
        switch brc
            case 'cau'
                divider = repmat(Parameters.LengthScales(CellIndex),[size(DistRC,1),1]);
                BiasRC(DistRC>0 & ioi) = (DistRC(DistRC>0 & ioi)-divider(DistRC>0 & ioi))/500;
                BiasRC(DistRC<0 & ioi) = (DistRC(DistRC<0 & ioi)-divider(DistRC<0 & ioi))/500;
            case 'ro'
                divider = repmat(Parameters.LengthScales(CellIndex),[size(DistRC,1),1]);
                BiasRC(DistRC>0 & ioi) = (DistRC(DistRC>0 & ioi)+divider(DistRC>0 & ioi))/500;
                BiasRC(DistRC<0 & ioi) = (DistRC(DistRC<0 & ioi)+divider(DistRC<0 & ioi))/500;
            case 'bi'
                BiasRC(:,ioi) = (abs(DistRC(:,ioi))-Parameters.LengthScales(CellIndex(ioi)))/500;
            case 'loc'
                divider = repmat(Parameters.LengthScales(CellIndex),[size(DistRC,1),1]);
                BiasRC(DistRC>0 & ioi) = (DistRC(DistRC>0 & ioi)./divider(DistRC>0 & ioi));
                BiasRC(DistRC<0 & ioi) = (DistRC(DistRC<0 & ioi)./divider(DistRC<0 & ioi));
        end
        BiasRC(DistRC==0 & ioi) = 0;
    end
    
    for bml = bML
        ioi = ismember(SC.Type,Parameters.Types(Parameters.BiasML == bml))';
        switch bml
            case 'lat'
                divider = repmat(Parameters.LengthScales(CellIndex),[size(DistML,1),1]);
                BiasML(DistML>0 & ioi) = abs(DistML(DistML>0 & ioi))./1000;
                BiasML(DistML<0 & ioi) = abs(DistML(DistML<0 & ioi))./100;
            case 'med'
                divider = repmat(Parameters.LengthScales(CellIndex),[size(DistML,1),1]);
                BiasML(DistML>0 & ioi) = abs(DistML(DistML>0 & ioi))./100;
                BiasML(DistML<0 & ioi) = abs(DistML(DistML<0 & ioi))./1000;
            case 'loc'
                BiasML(:,ioi) = abs(DistML(:,ioi))./1000;
        end
        BiasML(DistML==0 & ioi) = 0;
    end

    for bdv = bDV
        ioi = ismember(SC.Type,Parameters.Types(Parameters.BiasDV == bdv))';
        switch bdv
            case 'dor'
                divider = repmat(Parameters.LengthScales(CellIndex),[size(DistDV,1),1]);
                BiasDV(DistDV>0 & ioi) = abs(DistDV(DistDV>0 & ioi))./1000;
                BiasDV(DistDV<0 & ioi) = abs(DistDV(DistDV<0 & ioi))./100;
            case 'ven'
                divider = repmat(Parameters.LengthScales(CellIndex),[size(DistDV,1),1]);
                BiasDV(DistDV>0 & ioi) = abs(DistDV(DistDV>0 & ioi))./100;
                BiasDV(DistDV<0 & ioi) = abs(DistDV(DistDV<0 & ioi))./1000;
            case 'loc'
                BiasDV(:,ioi) = abs(DistDV(:,ioi))./1000;
        end
        BiasDV(DistDV==0 & ioi) = 0;
    end

    for bci = bCI
        ioi = ismember(SC.Type,Parameters.Types(Parameters.BiasContraIpsi == bci))';
        switch bci
            case 'contra'
                BiasCI(DistCI<0 & ioi) = 1;
                BiasCI(DistCI>0 & ioi) = 0;
            case 'ipsi'
                BiasCI(DistCI>0 & ioi) = 1;
                BiasCI(DistCI<0 & ioi) = 0;
            case 'bi'
                BiasCI(:,ioi) = 1;
        end
    end

    for bl = bLay
        if(~isundefined(bl{1}))
            ibl = cellfun(@(x) isempty(setxor(bl{1},x)),Parameters.BiasLayer);
            ioi = ismember(SC.Type,Parameters.Types(ibl))';
            BiasLay((~sum(contains(string(SC.Layers),string(bl{1})),2))&ioi) = 0.05;
        end
    end
    for bS = bSeg
        if(~isundefined(bS{1}))
            ibS = cellfun(@(x) isempty(setxor(bS{1},x)),Parameters.BiasSegment);
            ioi = ismember(SC.Type,Parameters.Types(ibS))';
            BiasSeg((~sum(contains(string(SC.Layers(:,4)),string(bS{1})),2))&ioi) = 0.5;
        end
    end

   for bmn = bMn
        if(~isundefined(bmn))
            imn = ismember(Parameters.BiasMN,bmn);
            ioi = ismember(SC.Type,Parameters.Types(imn))';
            if(contains(string(bmn),'only'))
                bmnn = erase(string(bmn),'only');
                BiasMn(((~ismember(SC.FlexExtID,bmnn))&ismember(SC.Type,'MN'))&ioi) = 0.05;
            else
                BiasMn(((~(ismember(SC.FlexExtID,bmn)|ismember(SC.FlexExtID,'Bi')))&ismember(SC.Type,'MN'))&ioi) = 0.3;
            end
        end
   end

   for bmn = bMn
        if(~isundefined(bmn))
            imn = ismember(Parameters.BiasMN,bmn);
            ioi = ismember(SC.Type,Parameters.Types(imn))';
            BiasModule(~(ioi')&ioi) = 0.75;
        end
    end

    Prob = zeros(size(BiasRC));
    if(all(isnan(BiasML(:))) && all(isnan(BiasDV(:))) && all(isnan(BiasRC(:))))
         Prob(:) = BiasCI(:).*BiasPop(:).*BiasLay(:).*BiasMn(:);
    else
         Prob(:) = geomean([normpdf(BiasML(:),0,1)./normpdf(0,0,1),normpdf(BiasDV(:),0,1)./normpdf(0,0,1),normpdf(BiasRC(:),0,1)./normpdf(0,0,1)],2,'omitnan').*BiasCI(:).*BiasPop(:).*BiasLay(:).*BiasSeg(:).*BiasMn(:);
    end
    Prob(isnan(Prob)) = 0;
    %% Add Random Component to connectivity if specified 
    if RC
        RandComp = zeros(size(Prob));
        ixs = randsample(1:1:numel(RandComp),round(numel(RandComp)/2),0);
        RandComp(ixs) = normrnd(0,RCvar,[1 numel(ixs)]);
        Prob = Prob + RandComp;
        Prob(Prob>1) = 1;
        Prob(Prob<0) = 0;
    end

    if SP
        Spar = zeros(size(Prob));
        ixs = randsample(1:1:numel(Spar),round(numel(Spar)*sparsity),0);
        Prob(ixs) = 0;
    end
    
    %% Compute Connection Probability and sparsify
    Sparsifier = binornd(1,Prob); % Sparsify based on binomial proba
    Connectivity = normrnd(1,var,size(Prob)); % Attribute random strength with mean 1
    Connectivity(:,E) = abs(Connectivity(:,E)).*Parameters.SynStrengths(CellIndex(E)); % Bias Mean 
    Connectivity(:,I) = -(Connectivity(:,I)).*Parameters.SynStrengths(CellIndex(I));
    
    Connectivity = Connectivity.*Sparsifier;
    Connectivity = Connectivity - diag(diag(Connectivity)); % Remove self synapse
    Connectivity = UnifySensoryInputProjection(Connectivity,SC);

    if(Bal)
        Connectivity(SC.Latera>0,:) = BalanceConnectivity(Connectivity(SC.Latera>0,:));
        Connectivity(SC.Latera<0,:) = BalanceConnectivity(Connectivity(SC.Latera<0,:));
    end
    if(resc)
        Connectivity(SC.Latera>0,:) = RescaleConnectivity(Connectivity(SC.Latera>0,:));
        Connectivity(SC.Latera<0,:) = RescaleConnectivity(Connectivity(SC.Latera<0,:));
    end
   
    obj.Sparsity = nnz(~Connectivity)/numel(Connectivity);
    obj.ConnMat = Connectivity;
    obj.Delays = DelaysPop;
    obj.Position = SC.Position;
    obj.Types = SC.Type;
    obj.Latera = SC.Latera;
    obj.Transmit = SC.Transmit;
    obj.Layers = categorical(SC.Layers(:,1));
    obj.Segment = categorical(SC.Layers(:,4));
    obj.MnID = SC.MnID;
    obj.FlexExtID = SC.FlexExtID;
    obj.ComputeEigenModesandNullSpace('Verbose',verbose)
end
%%
function Connectivity = UnifySensoryInputProjection(Connectivity,SC)
    sens = contains(string(SC.Type),'NF') | contains(string(SC.Type),'NP') | contains(string(SC.Type),'PEP') | contains(string(SC.Type),'TH');
    whr = logical(Connectivity);
    intype = repmat(SC.MnID',[length(SC.MnID) 1]);    
    intype(~whr) = '<undefined>';
    for jj = find(sens)'
        temp_out = intype(jj,:);
        temp_out(isundefined(temp_out)) = [];
        temp_out = unique(temp_out);
        if(~isempty(temp_out))
            mde = randsample(temp_out,1,true,ones(1,length(temp_out))/length(temp_out));
            %mde = mode(temp_out);
            notMN = (~(SC.MnID == mde))&(SC.Type=='MN');
            Connectivity(jj,notMN) = 0;
            Connectivity(notMN,jj) = 0;
        end
    end
end

%%



