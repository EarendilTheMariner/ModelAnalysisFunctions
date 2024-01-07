
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

%% Save vars in unique folder
for ii=1:100
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);
    [W, R] = BalancedRecRandom();
    close all
    clearvars folderName
    save('output.mat');
    cd ..    
end

for ii=1:100
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load("output.mat");
    PhaseDistance(W);
    close all
    clearvars folderName
    cd ..    
end
