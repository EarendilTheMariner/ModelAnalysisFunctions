h = figure;
for ii = 1:length(inhLoop)
        excLoop = linspace(LowerExcLS,inhLoop(ii),ExLSStep);
        for i = length(inhLoop)
            subplot(2,3,ii);
            scatter(excLoop,Summarys(ii,:));
            xlabel('V2a-2 (bidrectional bias) LS');
            ylabel('Oscillation amplitude (RMS)');
            legend(['V1 (bidirectional bias) LS = ',num2str(inhLoop(ii))]);
        end
end
savefig(h, 'Visualization.fig')