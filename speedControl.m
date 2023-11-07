function MouseNet = speedControl(varargin)

%% Oscillation frequency, input dependent recruitment of populations
    
    NP = NetworkParameters();
    NP.AddCelltype('V2a-1','cau','','','','','',1000,0.1000*5);
    NP.AddCelltype('V2a-2','bi','','','','','',500,0.1000*2.5);
    NP.AddCelltype('V1','bi','','','','','',2000,0.1000*2.5);
    MouseNet = Network(3000,0.5,NP);
    
    W = MouseNet.ConnMat;

    simTime = 30000;
    lowerIe = 20;
    upperIe = 30;
    V2a_1Thr = 25;
    V2a_2Thr = 22.5;
    V1_Thr = 20;

        for ii = 1:2:length(varargin)
            switch varargin{ii}
                case 'SimTime'
                    simTime = varargin{ii+1};
                case 'LowerInput'
                    lowerIe = varargin{ii+1};
                case 'UpperInput'
                    upperIe = varargin{ii+1};
                case 'V2a_1Threshold'
                    V2a_1Thr = varargin{ii+1};
                case 'V2a_2Threshold'
                    V2a_2Thr = varargin{ii+1};
                case 'V1Threshold'
                    V1_Thr = varargin{ii+1};
            end
         end

    V2a_1Ix = MouseNet.Types == 'V2a-1';
    V2a_2Ix = MouseNet.Types == 'V2a-2';
    V1Ix = MouseNet.Types == 'V1';

    threshold = ones(size(W,1),1);
    threshold(V2a_1Ix) = V2a_1Thr;
    threshold(V2a_2Ix) =  V2a_2Thr;
    threshold(V1Ix) = V1_Thr;

    inputSignal = zeros(size(W,1),simTime);
    inputSignal(:,1:5000) = lowerIe;
    inputSignal(:,5000:9999) = repmat(linspace(lowerIe,upperIe,5000),size(W,1),1);
    inputSignal(:,10000:20000) = upperIe;
    inputSignal(:,20000:24999) = repmat(linspace(upperIe,lowerIe,5000),size(W,1),1);
    inputSignal(:,25000:30000) = lowerIe;

    MouseNet.Simulate(simTime,'I_e', inputSignal, 'threshold', threshold);

    GeneratePlot(MouseNet);

end