
%% Simulation Heatmap
Rates = N.Rates;
[~,idx] = sort(N.Position(:,2));
Stats = zeros(6000,3088);
for i = 1:6000
    CurrRates = Rates(i,:);
    SortedRates = CurrRates(idx);
    Stats(i,:) = SortedRates;
end

figure;
imagesc(Stats);
set(gca, 'YDir', 'normal');
xlabel('Rostro-caudal position');
ylabel('Simulation time (ms)');
colorbar

%% Save vars in unique iteration folder
for ii=1:32
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii+68);
    mkdir(folderName);
    cd(folderName);
    [W, R] = BalancedRecRandom();
    close all
    clearvars folderName
    save('output.mat');
    cd ..    
end


%% Fix dynamics sim

for ii=1:68

    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load("output.mat");
    R = R(9001:end,:);
    clearvars folderName
    save("output.mat");
    delete Coords.mat
    delete PCVar.mat
    cd ..    
end


%% Compute kernel and save .png file in each iteration folder
for ii=1:100

    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load("output.mat");
    delete PCVar.mat
    delete kernel.png
    delete Coords.mat
    delete dynamics.png

    PhaseDistance(W);
    load('Rates.mat')
    % Eigenspectrum
    bevr = eig(W);
    [brevr,bixev] = sort(real(bevr),'descend');
    bievr = imag(bevr(bixev));

    % PCA
    [~, scores] = pca(Rates);
   
    % Plots 
    figure;
    subplot(2,3,2);
    scatter(brevr,bievr);
    vline(1);
    xlim([-1.5,1.5]);
    xlabel('Real');
    ylabel('Imaginary')
    title('Eigenspectrum');
    axis equal
    box off
    
    subplot(2,3,4:6);
    plot(Rates);
    grid();
    xlabel('Time (ms)');
    ylabel('Rate (Hz)');
    xlim([0,9999]);
    ylim([0,60]);
    box off

    subplot(2,3,1);
    imagesc(W);
    title('Connectivity matrix - Sparsity ≈ 0.9');
    box off

    subplot(2,3,3);
    scatter(scores(:,1),scores(:,2));
    xlabel('PC1');
    ylabel('PC2');
    xlim([-500,500]);
    ylim([-400,400]);
    box off

    set(gcf, 'WindowState', 'maximized');

    drawnow;  % Ensure the plot is fully rendered
    frame = getframe(gcf);  % Get the current frame
    im = frame2im(frame); % Convert the frame to an image
    imwrite(im, 'dynamics.png');


    close all
    clearvars folderName
    cd ..    

end


%% Compute number of oscillating networks 
NumOsc_N500 = 0;
for ii=1:100
    clearvars -except ii NumOsc_N500 NumOsc_N1000
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    if isfile('kernel.png')
        NumOsc_N500 = NumOsc_N500 + 1;
        clearvars folderName
    else
        clearvars folderName
    end
    cd ..    
end


