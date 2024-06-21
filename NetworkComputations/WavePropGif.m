function WavePropGif(Rates)

    for ii = 9000:10:14000
        close all
        clearvars -except Rates ii
        Space = Rates(ii,:);
        figure; 
        set(gcf, 'WindowState', 'maximized');
        bar(Space);
        ylim([0 70]);
        box off
    
        drawnow;  % Ensure the plot is fully rendered
        frame = getframe(gcf);  % Get the current frame
        im = frame2im(frame); % Convert the frame to an image
        imwrite(im, 'temp.png');
        im = imread('temp.png');
        [imind,cm] = rgb2ind(im,256);  % Convert the image to indexed color
    
        filename = 'WavePropagationC2.gif';
    
        % Write to GIF
        if ii == 9000
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.1);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.1);
        end
    end


end

