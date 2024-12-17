classdef LineNetwork < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Coordinates         % spatial coordinates of each neuron, length(coords) = network length, keeping neuron density fixed
        Rates               % time x n rate matrix
        Voltage             % time x n input potential matrix
        ConnMat             % n x n connectivity matrix
        E                   % n x 1 logical index of excitatory units
        I                   % n x 1 logical index of inhibitory units
        PDFs                % spatial projection pdfs for the five subpopulations 
        Projectome          % realized spatial projection profile of network over the pair-wise distance domain
        Size                % size (n) of network, equal to length to fix density
        Length              % spatial length of network
        PopParams           % 5 x 1 vector mean (Gaussian) projection lengths for each of the 5 subpopulations over the normalized pair-wise distance domain   
        Types               % n x 5 logical indices for subpopulations
        EigenMode           % 1st eigenvector of connectivity matrix
        Skewness            % Quantifies the skewness of the projectome. Defined as the median excitation distance. Distance where cumulative excitation reaches 0.5. 
        
     
    end
    
    methods
        function obj = LineNetwork(Length,PopParams)

            obj.PopParams = PopParams;
            obj.Length = Length;
            obj.Size = Length;
            obj.Coordinates = linspace(0,Length,Length);
            
            if length(PopParams)==5
                obj = Instantiate(obj); % Builds the network from basic parameters
            else
                obj = Instantiate2(obj);
            end

            obj.Projectome = SpatialCoupling(obj.ConnMat,obj.E,obj.I,obj.Coordinates,false); % Computes the spatial projectome of the network, 
     %       obj.Skewness = SkewnessFactor(obj.Projectome); 
            obj.Rates = [];      
            obj.EigenMode = EigenStructure(obj.ConnMat,1); % Computes leading eigenvector of connectivity matrix
        end

        
    end
    
    methods (Access = public)

        SimulateLine(obj,t_steps,varargin); % Simulates network, t_steps is duration of sim in milliseconds 
        GainMod(obj,Pop,Gain,normalizeFlag); % Scales synaptic weights of population of interest (Pop = 1 , ... , 5) by scaling facter Gain. 
    end

end

