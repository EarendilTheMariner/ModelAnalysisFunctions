function [SortingIndices,varargout] = GetFiringPhaseSorting(FiringMatrix,Ops,varargin)
    
    Method = 'Correlation';
    %Vars = var(FiringMatrix);
    %[~,I] = maxk(Vars,500);
    %Source = sum(FiringMatrix(:,I),2);
    [~, scores] = pca(FiringMatrix);
    Source = scores(:,1);
    for ii = 1:2:length(varargin)
        switch(varargin{ii})
            case 'Source'
                Source = varargin{ii+1};
            case 'Method'
                Method = varargin{ii+1};
        end
    end

    switch Method 
        case 'Correlation'
           % whole = FiringMatrix(:,15);
            whole = Source; 
            scores = [];
            corrs = [];
            for n = 1:size(FiringMatrix,2)
                [C,L] = xcorr(whole,FiringMatrix(:,n));
                [maxVal,loc] = max(C);
                scores = [scores L(loc)];
                corrs = [corrs maxVal];
            end
            
            [autocor,~] = xcorr(whole);
            [~,lcsh] = findpeaks(autocor);
            short = mean(diff(lcsh));
           

            scores =rem((scores.*2*pi)./short,2*pi);
            corrs = corrs./max(autocor);
            
            [scoresS,CorrSorting] =  sort(scores,'ascend');
            SortingIndices = CorrSorting;
            varargout = {scores,corrs};
            
        case 'Coherence'                
            %% Get Series with most rythmical activity
            Cohr = [];
            CohAmp = [];
            for ii = 1:size(FiringMatrix,2)
                [CohAm,CohA,~] = ComputeCoherence(Source,FiringMatrix(:,ii),10,1/Ops.fs,5);
                Cohr = [Cohr CohA];
                CohAmp = [CohAmp CohAm];
            end
            scores =rem((Cohr.*2*pi)./range(Cohr),2*pi);
            [Cohrs,CohSorting] =  sort(Cohr,'ascend');
            SortingIndices = CohSorting;
            varargout = {Cohr,CohAmp};
    end
end
