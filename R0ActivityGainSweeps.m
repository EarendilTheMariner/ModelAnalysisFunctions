function Stats = R0ActivityGainSweeps(varargin)

    ActGainSweep = linspace(0,2,5);
    Summarys = zeros(1,length(ActGainSweep));
    

    NetParams = NetworkParameters();
    NetParams.AddCelltype('V2a-1','cau','','','','','',1000,0.1*1.8);
    NetParams.AddCelltype('V1','bi','','','','','',3300,0.1*1.8);
    MouseNet = Network(3000,0.5,NetParams);
    W = MouseNet.ConnMat;
         
    for i = 1:length(ActGainSweep)
        close all 
        clearvars -except Summarys i ActGainSweep MouseNet W
   
        MouseNet.Simulate(10000, 'gain', ones(size(W,1),10000).*ActGainSweep(i));

        [~, scores] = pca(MouseNet.Rates);
        pc1Scores = scores(:,1);
        amplitude = rms(pc1Scores); % "Amplitude" from projected rates onto PC1 and taking root mean square (RMS)

        uniqueFilename = ['MouseNet_',num2str(i), '.mat']; 

        save(uniqueFilename, 'MouseNet'); %Saves the network in your current directory
        
        Summarys(1:i) = amplitude;
        
       
        figure;
        subplot(2,3,4:6);
        plot(MouseNet.Rates);
        xlabel('Time (ms)');
        ylabel('Firing rate');
        grid();

        subplot(2,3,1);
        scatter(ActGainSweep,Summarys)
        xlabel('Gain (Iteration)');
        ylabel('Amplitude (RMS)');
        grid();

        subplot(2,3,3);
        scatter(scores(:,1),scores(:,2))
        xlabel('PC1')
        ylabel('PC2')
        grid();

                
        drawnow;  % Ensure the plot is fully rendered
        frame = getframe(gcf);  % Get the current frame
        im = frame2im(frame);  % Convert the frame to an image
        [imind,cm] = rgb2ind(im,256);  % Convert the image to indexed color

        filename = 'ActGainSweep.gif'

        % Write to GIF
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.1);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.1);
        end



    end             
    
    save('Stats',"Summarys");

end