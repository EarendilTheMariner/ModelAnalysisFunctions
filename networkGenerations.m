function Summarys = networkGenerations(varargin)

%% Simulates networks and computes summary statistics 

    %Default upper and lower limits of excitatory and inhibitory length scales and
    %iteration resolution

    UpperExcLS = 5000;
    LowerExcLS = 300;
    ExLSStep = 100;

    UpperInhLS = 10000;
    LowerInhLS = 5000;
    InLSStep = 5;
    
    %variable arguments

    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'UpperLengthScaleEx'
                UpperExcLS = varargin{ii+1};
            case 'LowerLengthScaleEx'
                LowerExcLS = varargin{ii+1};
            case 'ExLengthScaleStep'
                ExLSStep = varargin{ii+1};
            case 'UpperLengthScaleInh'
                UpperInhLS = varargin{ii+1};
            case 'LowerLengthScaleInh'
                LowerInhLS = varargin{ii+1};
            case 'InhLengthScaleStep'
                InLSStep = varargin{ii+1};
        end
    end
    
    %Inner (excitatory) and outer (inhibitory) loop

    excLoop = linspace(LowerExcLS,UpperExcLS,ExLSStep);
    inhLoop = linspace(LowerInhLS,UpperInhLS,InLSStep);
    
    %Dimensions for summary stats (I've chosen five outer loop iterations)
    Summarys = zeros(3,length(excLoop));
    Summarys(:,:,2) = zeros(3,length(excLoop));
    Summarys(:,:,3) = zeros(3,length(excLoop));
    Summarys(:,:,4) = zeros(3,length(excLoop));
    Summarys(:,:,5) = zeros(3,length(excLoop));

    %Main simulation loop. For each inhibitory length scale, iterate
    %through exictatory length scales "LowerExcLS" to the current
    %inhibitory length scale. 

    for ii = 1:length(inhLoop)
        excLoop = linspace(LowerExcLS,inhLoop(ii),ExLSStep);
        for i = 1:length(excLoop)
            close all 
            clearvars -except Summarys excLoop inhLoop i ii LowerExcLS ExLSStep
       
            NetParams = NetworkParameters();
            NetParams.AddCelltype('V2a-1','cau','','','','','',excLoop(i),0.3850);
            NetParams.AddCelltype('V2b','bi','','','','','',inhLoop(ii),0.11);
            MouseNet = Network(3000,0.5,NetParams);
            
            MouseNet.Simulate(10000);
            
            GeneratePlot(MouseNet); %Plots gain-weighted projection distribution, eigenspectrum, simulated rates
            
            eigs = eig(MouseNet.ConnMat);
            maxEig = max(real(eigs));
            eigIx = real(eigs) == maxEig;
            eigFreq = imag(eigs(eigIx))/2*pi; %predicted frequency from dominant oscillatory mode
        
            LR = MouseNet.Parameters.LengthScales(1)/MouseNet.Parameters.LengthScales(2); % exc. to inh. lengthscale ratio
            
        
            [~, scores] = pca(MouseNet.Rates);
            pc1Scores = scores(:,1);
            amplitude = rms(pc1Scores); % "Amplitude" from projection rates onto PC1 and taking root mean square (RMS)
            col = [eigFreq(1), LR, amplitude];
        
            Summarys(:,i,ii) = col'; % Stores statistics (frequency, amplitude, LR) from current network

            uniqueFilename = ['MouseNet_', num2str(ii),'_',num2str(i), '.mat']; 

            save(uniqueFilename, 'MouseNet'); %Saves the network in your current directory


        end
        
    end

    % saves the stats from all the networks as 'Stats.mat' file
    save('Stats',"Summarys");

        h = figure;    
    for i = 1:5
        subplot(2,3,i);
        scatter(Summarys(2,:,i),Summarys(1,:,i));
        xlabel('Lengthscale ratio (Exc/Inh)');
        ylabel('Frequency (Hz)');
        legend(['Inhibitory length scale = ',num2str(((i-1)+5)*1000)]);
    end
    savefig(h, 'Visualization.fig')

end



