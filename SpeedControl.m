function MouseNet = SpeedControl(varargin)

%% Oscillation frequency, input dependent recruitment of populations
    
    NP = NetworkParameters();
    NP.AddCelltype('V2a-1','cau','','','','','',1000,0.1000*2.5);
    NP.AddCelltype('V2a-2','bi','','','','','',500,0.1000*2.5);
    NP.AddCelltype('V1','bi','','','','','',2000,0.1000*2.5);
    MouseNet = Network(3000,0.5,NP);
    
    W = MouseNet.ConnMat;

    simTime = 34000;
    lowerIe = 20;
    upperIe = 20;
    V2a_1Thr = 20;
    V2a_2Thr = 20;
    V1_Thr = 20;
    caudInp = false;

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
                case 'CaudalInput'
                    caudInp = varargin{ii+1};
            end
         end

    V2a_1Ix = MouseNet.Types == 'V2a-1';
    V2a_2Ix = MouseNet.Types == 'V2a-2';
    V1Ix = MouseNet.Types == 'V1';

    threshold = ones(size(W,1),1);
    threshold(V2a_1Ix) = V2a_1Thr;
    threshold(V2a_2Ix) = V2a_2Thr;
    threshold(V1Ix) = V1_Thr;

    inputSignal = ones(size(W,1),simTime)*lowerIe;

    if caudInp

        inputSignal(V2a_1Ix,1:3*simTime/10) = lowerIe;
        inputSignal(V2a_1Ix,3*simTime/10:(3.5*simTime/10)-1) = repmat(linspace(lowerIe,upperIe,simTime/20),sum(V2a_1Ix),1);
        inputSignal(V2a_1Ix,3.5*simTime/10:6.5*simTime/10) = upperIe;
        inputSignal(V2a_1Ix,6.5*simTime/10:(7*simTime/10)-1) = repmat(linspace(upperIe,lowerIe,simTime/20),sum(V2a_1Ix),1);
        inputSignal(V2a_1Ix,7*simTime/10:simTime) = lowerIe;
    
    else 

        inputSignal(:,1:3*simTime/10) = lowerIe;
        inputSignal(:,3*simTime/10:(3.5*simTime/10)-1) = repmat(linspace(lowerIe,upperIe,simTime/20),size(W,1),1);
        inputSignal(:,3.5*simTime/10:6.5*simTime/10) = upperIe;
        inputSignal(:,6.5*simTime/10:(7*simTime/10)-1) = repmat(linspace(upperIe,lowerIe,simTime/20),size(W,1),1);
        inputSignal(:,7*simTime/10:simTime) = lowerIe;
    end


    MouseNet.Simulate(simTime,'I_e', inputSignal, 'threshold', threshold);
    MouseNet.PlotRates('InputSignal', inputSignal);


end