function MouseNet = SpeedControl(varargin)

%% Oscillation frequency, input dependent recruitment of populations
  
    lowerIe = 10;
    upperIe = 25;
    V2a_1Thr = 25;
    V2a_2Thr = 20;
    V1_Thr = 20;

        for ii = 1:2:length(varargin)
            switch varargin{ii}
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


    CauThrSweep = [linspace(0,50,100)];
    NP = NetworkParameters();
    NP.AddCelltype('V2a-1','cau','','','','','',1000,0.5);
    NP.AddCelltype('V2a-2','bi','','','','','',500,0.18);
    NP.AddCelltype('V1','bi','','','','','',3300,0.18);
    MouseNet = Network(3300,0.5,NP);
    W = MouseNet.ConnMat;
    threshold = ones(size(W,1),1)*20;
    V2a_1Ix = MouseNet.Types == 'V2a-1';
    close all

    
    for i = 1:length(CauThrSweep)
        clearvars -except CauThrSweep i MouseNet V2a_1Ix W threshold
        close all
        MouseNet.Rates=[];
        
        threshold(V2a_1Ix) = CauThrSweep(i);
        MouseNet.Simulate(8000,'threshold', threshold);
        
        uniqueFilename = ['MouseNet_',num2str(i), '.mat']; 
        save(uniqueFilename, 'MouseNet');

        [~, scores] = pca(MouseNet.Rates);

        figure;
        subplot(2,2,3:4);
        plot(MouseNet.Rates);
        xlabel('Time (ms)');
        ylabel('Firing rate');
        xlim([0,8000])
        ylim([0, 50])
        grid();
        box off
    
        subplot(2,2,2);
        scatter(scores(:,1),scores(:,2))
        xlabel('PC1');
        ylabel('PC2');
        xlim([-500,500]);
        ylim([-400,400]);
        box off

        subplot(2,2,1);
        scatter(i,CauThrSweep(i));
        xlim([1,100]);
        ylim([0,50]);
        xlabel('Iteration');
        ylabel('Caudal (V2a-1) Threshold');
        grid();
        box off


        set(gcf, 'WindowState', 'maximized');
    
        drawnow;  % Ensure the plot is fully rendered
        frame = getframe(gcf);  % Get the current frame
        im = frame2im(frame); % Convert the frame to an image
        imwrite(im, 'temp.png');
        im = imread('temp.png');
        [imind,cm] = rgb2ind(im,256);  % Convert the image to indexed color
    
        filename = 'R1CauActGainSweep.gif';
    
        % Write to GIF
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.4);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.4);
        end
    end
    delete('temp.png');
end