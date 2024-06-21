classdef Network < handle
    properties
        ConnMat
        Sparsity 
        Delays
        Position 
        Types 
        Latera
        Transmit
        Layers  
        Segment
        MnID
        FlexExtID 
        Rates
        EstimatedRates
        Spikes
        PC
        Phase 
        NullSpace
        EigenValues
        EigenModes
        Parameters
    end

    properties (Hidden)
        Geometry
        Genetics 
    end

    methods (Access = public)
        function obj = Network(Length,Density,Symmetry,Parameters,varargin)
            obj.Parameters = Parameters;
            obj.Geometry = SpinalCordGeometry.GetMouseGeometry(Length,Density,Symmetry); 
            obj.Genetics = SpinalCordGenetics.GetMouseGenetics();
            obj.InstantiateNetwork(obj.Geometry,obj.Parameters,varargin{:});
        end
         
    end

    methods (Access = public)
        Simulate(obj,t_steps,varargin)
    end

    methods(Access=public)
        InstantiateNetwork(obj,Geometry,Parameters,varargin)
        PlotRates(obj,varargin)
        PlotRaster(obj,varargin)
        PlotEMG(obj,varargin)
        PlotNeurogram(obj,varargin)
        AnimatedPlotRates(obj,varargin)
        ComputeEigenModesandNullSpace(obj,varargin)
    end

    methods(Access=public)
        function  ComputePhase(obj,varargin)
           
            UoI = true(size(obj.ConnMat,1),1);
            for ii = 1:2:length(varargin)
                switch varargin{ii}
                    case 'UoI'
                       UoI = varargin{ii+1};
                end
            end
            FiringMat = obj.Rates(:,UoI);
            [~,s,~] = pca(FiringMat);
            whole = s(:,1); 
            scores = [];
            corrs = [];
            for n = 1:size(obj.Rates,2)
                [C,L] = xcorr(whole,obj.Rates(:,n));
                [maxVal,loc] = max(C);
                scores = [scores L(loc)];
                corrs = [corrs maxVal];
            end
            
            [autocor,~] = xcorr(whole);
            [~,lcsh] = findpeaks(autocor);
            short = mean(diff(lcsh));
           
            scores =rem((scores.*2*pi)./short,2*pi);
            obj.Phase = scores';
        end
        function  ComputePC(obj,varargin)
            UoI = true(size(obj.ConnMat,1),1);
            Est = 0;
            Source = true(1,size(obj.Rates,1));
            for ii = 1:2:length(varargin)
                switch varargin{ii}
                    case 'UoI'
                       UoI = varargin{ii+1};
                    case 'Estimated'
                       Est = varargin{ii+1};
                    case 'Source'
                        Source = varargin{ii+1};
                end
            end
            if(Est)
                FiringMat = obj.EstimatedRates;
            else
                FiringMat = obj.Rates;
            end
            FM = FiringMat;
            [c,s,~] = pca(FM(Source,UoI));
            obj.PC = FiringMat(:,UoI)*c;
        end
        function  ComputeEstimatedRates(obj,varargin)
            if(isempty(obj.Spikes))
                obj.ComputeSpikes;
            end
            obj.EstimatedRates = GetGaussianFiring(obj.Spikes,50,1000); 
        end
        function  ComputeSpikes(obj,varargin)
            UoI = true(size(obj.ConnMat,1),1);
            for ii = 1:2:length(varargin)
                switch varargin{ii}
                    case 'UoI'
                       UoI = varargin{ii+1};
                end
            end
            RI =  obj.Rates(1:end,UoI);
            RIsp = ((RI-min(RI,[],1)));
            RIsp = RIsp./(max(RIsp,[],1));
            Poiss = poissrnd(RIsp/50,size(RIsp));
            obj.Spikes = logical(Poiss);
        end

    end
end