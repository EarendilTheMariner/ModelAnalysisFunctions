function Coords = PhaseDistance(W,varargin)
    
        n = length(W(:,1));
        f_exc = 0.5;

        for ii = 1:2:length(varargin)
            switch varargin{ii}
                case 'f_E'
                    f_exc = varargin{ii+1};
            end
        end

        % Extract representative signal with PCA
        Rates = SimulateNetwork(W,20000);
        Rates = Rates(10001:end-1,:);                                           
        [~,scores] = pca(Rates);
        PCRef = scores(:,1);

        % Zero padding to reduce boundary effects
        PadLength = 10000;
        StartPad = zeros(PadLength,1);
        EndPad = StartPad;

        % Mirror padding
        StartPad = flipud(StartPad);
        EndPad = flipud(EndPad);
        
        % Padded PCRef      
        PCRef = [StartPad; PCRef; EndPad];

        % Smooth padded signal and trim pads

        PCRef = smoothdata(PCRef);
        PCRef = PCRef(PadLength+1:end-PadLength);

        %Identify peaks and compute variance

        [~,locs] = findpeaks(PCRef); 
        PCVar = var(PCRef);

        if PCVar > 1000 & length(locs) > 1

            Per = mean(diff(locs));
            
            PeriodRates = Rates(length(Rates(:,1))*0.5-0.5*Per:length(Rates(:,1))*0.5+0.5*Per,:);
    
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