loop = linspace(0.1,1.5,100);

for i = 1:length(loop)   
        clearvars MouseNet
        load(['MouseNet_',num2str(i),'.mat'])
        [~, scores] = pca(MouseNet.Rates);
        pc1Scores = scores(:,1);
        amplitude = rms(pc1Scores); % "Amplitude" from projected rates onto PC1 and taking root mean square (RMS)

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
        scatter(loop(i)*1000,Summarys(1:i))
        xlim([0,2000]);
        ylim([0,500]);
        xlabel('Caudal Length Scale (Iteration)');
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

        filename = 'CauLSSweepHD.gif';

        % Write to GIF
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.4);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.4);
        end
        close all

end             
delete('temp.png');