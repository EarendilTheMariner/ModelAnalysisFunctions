Rates = MouseNet.Rates;
[~,idx] = sort(MouseNet.Position(:,2));
Stats = zeros(10000,689);
for i = 1:10000
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

