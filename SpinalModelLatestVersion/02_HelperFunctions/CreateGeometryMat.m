            SN = categorical({'TH','PEP1','PEP2','NP1','NP2','NP3','NF1','NF2','NF3','NF4','NF5'})';
            % Load necessary data
            load('01_SpinalGeometryFiles\matResults.mat')
            load('01_SpinalGeometryFiles\matCoord.mat')
            load('01_SpinalGeometryFiles\matGeneCellTypes.mat')
            load("01_SpinalGeometryFiles\matLayers.mat")
            load("01_SpinalGeometryFiles\matSegmentDimensions.mat")

            n = 1000;
            rightpos = [max(RCTDCoord.x)+range(RCTDCoord.x)/5,min(RCTDCoord.y)+range(RCTDCoord.y)/10];
            leftpos = [min(RCTDCoord.x)-range(RCTDCoord.x)/5,min(RCTDCoord.y)+range(RCTDCoord.y)/10];

            for jj = 1:length(SN)
                RCTDCellTypes = cat(1,RCTDCellTypes,{'',string(SN(jj)),1});
                weights(size(weights,1)+[1:(n/2)],size(weights,2)+1) = randsample(maxk(unique(weights(:)),1000),n/2);
                weights(size(weights,1)+[1:(n/2)],size(weights,2)) = randsample(maxk(unique(weights(:)),1000),n/2,1);
                matLayers = [matLayers;repmat({'DRG','DRG','DRG',rightpos(1),rightpos(2)},[n/2,1])];
                matLayers = [matLayers;repmat({'DRG','DRG','DRG',leftpos(1),leftpos(2)},[n/2,1])];
                RCTDCoord = [RCTDCoord;repmat({'',rightpos(1)+normrnd(0,10,1),rightpos(2)},[n/2,1])];
                RCTDCoord = [RCTDCoord;repmat({'',leftpos(1)+normrnd(0,10,1),leftpos(2)},[n/2,1])];                
            end
            matLayers(end-n+1:end,[4 5]) =  matLayers(end-n+1:end,[4 5]) + normrnd(0,10,[n 2]);
            RCTDCoord(end-n+1:end,[2 3]) =  RCTDCoord(end-n+1:end,[2 3]) + normrnd(0,10,[n 2]);


            clear jj leftpos results_df rightpos score_mat SN weights_doublet;
            save("03_Utilities\Spinal Geometry\matGeometry.mat","RCTDCoord","matLayers","RCTDCellTypes","Segments","SegmentLength","weights",'-mat');
            clear