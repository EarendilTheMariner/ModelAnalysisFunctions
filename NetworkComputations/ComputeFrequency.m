function [Per, Frequency] = ComputeFrequency(Rates)

    %% 

    
    if isempty(Rates)
        Per = NaN;
        Frequency = NaN;
    else
        [~,scores] = pca(Rates);
        PCRef = scores(:,1);
        [~,locs] = findpeaks(PCRef,'MinPeakProminence',5);
        
        if length(locs) > 1 & var(PCRef) > 1
    
        Per = mean(diff(locs));
        Frequency = 1./(Per/1000);

        else    
        
        Per = NaN;
        Frequency = NaN;
        end
         
    end

end


     


