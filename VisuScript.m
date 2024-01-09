
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

%% Compute kernel and save .png file in each iteration folder
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

for i = 1:500
    for j = 1:500
        SymDist(i,j) = SymDist(i,j)*1000;
    end
end
