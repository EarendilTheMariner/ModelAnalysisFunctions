function Stats = R0SynGainSweep(varargin)

    SynGainSweep = linspace(0.5,5,50);
    Summarys = zeros(1,length(SynGainSweep));
    Caudal = false;
    Bidi = false;

    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'Caudal'
                Caudal = varargin{ii+1};
            case 'Bidirectional'
                Bidi = varargin{ii+1};
        end
    end
    
         
    for i = 1:length(SynGainSweep)
        close all 
        clearvars -except Summarys i SynGainSweep MouseNet W Caudal Bidi

        if Caudal

            NetParams = NetworkParameters();
            NetParams.AddCelltype('V2a-1','cau','','','','','',1000,0.1*SynGainSweep(i));
            NetParams.AddCelltype('V1','bi','','','','','',3300,0.1*1.8);
            MouseNet = Network(3300,0.5,NetParams);   
            MouseNet.Simulate(10000);

        elseif Bidi

            NetParams = NetworkParameters();
            NetParams.AddCelltype('V2a-1','cau','','','','','',1000,0.1*1.8);
            NetParams.AddCelltype('V1','bi','','','','','',3300,0.1*SynGainSweep(i));
            MouseNet = Network(3000,0.5,NetParams);   
            MouseNet.Simulate(10000);

        else
            NetParams = NetworkParameters();
            NetParams.AddCelltype('V2a-1','cau','','','','','',1000,0.1*SynGainSweep(i));
            NetParams.AddCelltype('V1','bi','','','','','',3300,0.1*SynGainSweep(i));
            MouseNet = Network(3000,0.5,NetParams);   
            MouseNet.Simulate(10000);
        end

    
        [~, scores] = pca(MouseNet.Rates);
        pc1Scores = scores(:,1);
        amplitude = rms(pc1Scores); % "Amplitude" from projected rates onto PC1 and taking root mean square (RMS)

        uniqueFilename = ['MouseNet_',num2str(i), '.mat']; 

        save(uniqueFilename, 'MouseNet'); %Saves the network in your current directory
        
        Summarys(1:i) = amplitude;
        
        close all
        
       
        figure;
        subplot(2,2,3:4);
        plot(MouseNet.Rates);
        xlabel('Time (ms)');
        ylabel('Firing rate');
        xlim([0,10000])
        ylim([0, 50])
        grid();
        box off

        subplot(2,2,1);
        scatter(SynGainSweep(i)*0.1,Summarys(1:i))
        xlim([0,5*0.1]);
        ylim([0,500]);
        xlabel('Caudal Synaptic Gain (Iteration)');
        ylabel('Amplitude (RMS)');
        grid();
        box off

        subplot(2,2,2);
        scatter(scores(:,1),scores(:,2))
        xlabel('PC1');
        ylabel('PC2');
        xlim([-500,500]);
        ylim([-400,400]);
        box off

                
        set(gcf, 'WindowState', 'maximized');

        drawnow;  % Ensure the plot is fully rendered
        frame = getframe(gcf);  % Get the current frame
        im = frame2im(frame); % Convert the frame to an image
        imwrite(im, 'temp.png');
        im = imread('temp.png');
        [imind,cm] = rgb2ind(im,256);  % Convert the image to indexed color

        filename = 'CauSynGainSweep.gif'

        % Write to GIF
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.4);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.4);
        end
    end

    delete('temp.png');
    save('Stats',"Summarys");
end