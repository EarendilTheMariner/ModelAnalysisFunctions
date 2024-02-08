function Coords = PhaseDistance(W,Rates,varargin)
    
    n = length(W(:,1));

    % Extract representative signal with PCA

    [~,scores] = pca(Rates);
    PCRef = scores(:,1);

    %Identify peaks and compute variance

    [~,locs] = findpeaks(PCRef,'MinPeakProminence',5);

    if length(locs) > 1

        Per = mean(diff(locs));
        
        PeriodRates = Rates(1+length(Rates(:,1))*0.5-0.5*Per:length(Rates(:,1))*0.5+0.5*Per,:);

        NPeaks = zeros(n,1);

        for i = 1:length(W(:,1))
            [~,loc] = max(PeriodRates(:,i));
            NPeaks(i) = loc;
        end

        Coords = (NPeaks-min(NPeaks))/max(NPeaks);
        Coords = Coords*1000;  

    else
        Coords = [];
        disp('No or too slow oscillations')
    end   
end