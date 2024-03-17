function WavePropGif(Rates,W,I,E,Coords)

    [Per,~] = ComputeFrequency(Rates);
    PerRates = Rates(3000:3000+Per,:);
    Snaps = round(linspace(1,Per,100));
    
    for ii = 1:length(Snaps)
        close all
        clearvars -except Snaps ii PerRates Per W E I Coords
        Space = PerRates(Snaps(ii),:);
        figure; 
        subplot(1,2,2)
        set(gcf, 'WindowState', 'maximized');
        bar(Space);
        ylim([0 45]);
        xlabel('1D Space','FontSize',20);
        ylabel('Firing Rate', 'FontSize',20);
        subplot(1,2,1)
        SpatialCoupling(W,E,I,Coords);
    
    
        drawnow;  % Ensure the plot is fully rendered
        frame = getframe(gcf);  % Get the current frame
        im = frame2im(frame); % Convert the frame to an image
        imwrite(im, 'temp.png');
        im = imread('temp.png');
        [imind,cm] = rgb2ind(im,256);  % Convert the image to indexed color
    
        filename = 'WavePropagationC2.gif';
    
        % Write to GIF
        if ii == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.2);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.2);
        end
    end


end

