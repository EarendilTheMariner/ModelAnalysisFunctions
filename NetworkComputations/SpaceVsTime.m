function SpaceVsTime(Rates)

         
        Rates = downsample(Rates,50);

        fig = figure;
        imagesc(flipud(Rates));
        box off
        set(gca, 'YDir', 'normal'); 
        set(gcf, 'WindowState', 'maximized');
        colormap(gray);

        drawnow;  % Ensure the plot is fully rendered
        save2pdf(fig,['./'],['SpaceVsTime'],'-dpdf');


end
        


    