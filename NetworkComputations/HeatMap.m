function HeatMap(Rates,saveFlag)

    fig = figure; 
    set(gcf, 'WindowState', 'maximized');

    imagesc(downsample(Rates,10)');
    set(gca, 'YDir', 'normal');
    box off
    colormap(gray);
    colorbar
    drawnow;  

    xlabel('Time', 'FontSize',20);
    ylabel('Space','FontSize',20);
   

    if saveFlag
    save2pdf(fig,['./'],['HeatMap_Unscrambled_Global'],'-dpdf');
    end


end

