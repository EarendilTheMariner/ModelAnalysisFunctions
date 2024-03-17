function Stats = HatMotif(varargin)
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

    Summarys = zeros(length(inhLoop),length(excLoop));

    %Main simulation loop. For each inhibitory length scale, iterate
    %through exictatory length scales "LowerExcLS" to the current
    %inhibitory length scale. 

    for ii = 1:length(inhLoop)
        excLoop = linspace(LowerExcLS,inhLoop(ii),ExLSStep);
        for i = 1:length(excLoop)
            close all 
            clearvars -except Summarys excLoop inhLoop i ii LowerExcLS ExLSStep
       
            NetParams = NetworkParameters();
            NetParams.AddCelltype('V2a-1','cau','','','','','',1500,0.1*2.5);
            NetParams.AddCelltype('V2a-2','bi','','','','','',excLoop(i),0.1*2.5)
            NetParams.AddCelltype('V1','bi','','','','','',inhLoop(ii),0.1*2.5);
            MouseNet = Network(3000,0.5,NetParams);
            
            MouseNet.Simulate(10000);
            
            GeneratePlot(MouseNet); %Plots gain-weighted projection distribution, eigenspectrum, simulated rates
            
        
            [~, scores] = pca(MouseNet.Rates);
            pc1Scores = scores(:,1);
            amplitude = rms(pc1Scores); % "Amplitude" from projected rates onto PC1 and taking root mean square (RMS)


            uniqueFilename = ['MouseNet_', num2str(ii),'_',num2str(i), '.mat']; 

            save(uniqueFilename, 'MouseNet'); %Saves the network in your current directory
            
            Summarys(ii,i) = amplitude;
        end
        
    end

    % saves the stats from all the networks as 'Stats.mat' file
    save('Stats',"Summarys");

end


    