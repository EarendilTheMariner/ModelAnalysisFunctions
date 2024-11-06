for ii = 1:100
    close all
    clearvars -except ii 
 %   folderName = sprintf('Iteration_%d', ii);
 %   cd(folderName);    
 %   load("DiffCtrl.mat");
 %   load("EigenvaluesCtrl.mat");
    load(sprintf('Sp2_GainMod_%d', ii));
  %  load("Rates.mat");

%    cd ..

    
    figure;
    
    fig = gcf;
    
    subplot(2,2,1);
    SpatialCoupling(N.ConnMat,N.E,N.I,N.Coordinates,true);
    ylim([-50 50]);
    box off

    subplot(2,2,2);
    EigenSpectrum(N.ConnMat);
    axis equal
    box off

    
    subplot(2,2,3:4);
    plot(N.Rates(10001:end,:));
    ylim([0 70]);
    ylabel('Firing Rate', 'FontSize',20);
    xlabel('Time','FontSize',20);
    box off

    set(fig, 'WindowState', 'maximized');
    
    % Save the figure directly to a file
    print(fig, 'temp.png', '-dpng', '-r300');  % Save at 300 dpi
    
    im = imread('temp.png');
    [imind, cm] = rgb2ind(im, 256);  % Convert the image to indexed color
    
    filename = 'Overview.gif';
    
    % Write to GIF
    if ii == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end