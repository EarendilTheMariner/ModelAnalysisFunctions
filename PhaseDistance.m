function Coords = PhaseDistance(W,varargin)
    
        n = length(W(:,1));
        f_exc = 0.5;

        for ii = 1:2:length(varargin)
            switch varargin{ii}
                case 'f_E'
                    f_exc = varargin{ii+1};
            end
        end

        
        Rates = SimulateNetwork(W,20000);
        Rates = Rates(10000:end-1,:);                                           %Trim off transient initial dynamics
        SmoothR = smoothdata(Rates);                                            %Smooth high frequency fluctuation
        %[~,locs] = findpeaks(SmoothR(:,1));                                    %Peak times of arbitrary representative unit 
        [~,scores] = pca(SmoothR);
        PCRef = scores(:,1);
        [~,locs] = findpeaks(PCRef); 
        PCVar = var(PCRef);

        if PCVar > 1000 & length(locs) > 1

            Per = mean(diff(locs));
           %MidPeaks = [locs(round(length(locs)/2)),locs(round(length(locs)/2)-1)]; %Central peak-times of rep. unit
            
            
           %PeriodRates = SmoothR(MidPeaks(2):MidPeaks(1),:);
            PeriodRates = SmoothR(length(SmoothR(:,1))*0.5-0.5*Per:length(SmoothR(:,1))*0.5+0.5*Per,:);
    
            NPeaks = zeros(n,1);
    
            for i = 1:length(W(:,1))
                [~,loc] = max(PeriodRates(:,i));
                NPeaks(i) = loc;
            end
    
            XCoords = (NPeaks-min(NPeaks))/max(NPeaks);
            XCoords = XCoords*1000;
    
            fE = f_exc*n;
            fI = (1-f_exc)*n;
            E = false(1,n);
            I = E;
            E(1:fE) = true;
            I(fE+1:end) = true;
    
            Coords = XCoords;
            
            save('Coords.mat','Coords');
            
            CouplingKernel(W,f_exc,Coords,Rates);
            disp(var(PCRef))
            save('PCVar','PCVar');

        else
            disp(var(PCRef));
            disp('No oscillations');
            save('PCVar','PCVar');
        end
        
end