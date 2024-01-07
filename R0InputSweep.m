function Stats = R0InputSweep(varargin)

    InputSweep = linspace(0.1,40,80);
    Summarys = zeros(1,length(InputSweep));
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

    NetParams = NetworkParameters();
    NetParams.AddCelltype('V2a-1','bi','','','','','',1000,0.1*1.8);
    NetParams.AddCelltype('V1','bi','','','','','',3300,0.1*1.8);
    MouseNet = Network(3300,1,NetParams);
    W = MouseNet.ConnMat;

    V2a_1Ix = MouseNet.Types == 'V2a-1';
    V1Ix = MouseNet.Types == 'V1';
    BaseInput = ones(size(W,1),10000).*20;

         
    for i = 1:length(InputSweep)
        close all 
        clearvars -except Summarys i InputSweep MouseNet W V2a_1Ix V1Ix Caudal Bidi BaseInput
        
        MouseNet.Rates = [];

        if Caudal
            BaseInput(V2a_1Ix,:) = InputSweep(i);
            MouseNet.Simulate(10000, 'I_e', BaseInput);

        elseif Bidi 
            BaseInput(V1Ix,:) = InputSweep(i);
            MouseNet.Simulate(10000, 'I_e', BaseInput);

        else
            BaseInput(:,:) = InputSweep(i);
            MouseNet.Simulate(10000, 'I_e', BaseInput);
        end

        [~, scores] = pca(MouseNet.Rates);
        pc1Scores = scores(:,1);
        amplitude = rms(pc1Scores); % "Amplitude" from projected rates onto PC1 and taking root mean square (RMS)

        uniqueFilename = ['MouseNet_',num2str(i), '.mat']; 

        save(uniqueFilename, 'MouseNet'); %Saves the network in your current directory
        
        Summarys(1:i) = amplitude;
       
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
        scatter(InputSweep(i),Summarys(1:i))
        xlim([0,40]);
        ylim([0,500]);
        xlabel('Bidirectional (V1) Input (Iteration)');
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

        filename = 'BidiInputSweep.gif';

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



