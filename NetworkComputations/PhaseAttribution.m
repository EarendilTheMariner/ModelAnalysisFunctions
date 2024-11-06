function [Diffs, MeanProjLength] = PhaseAttribution(W,Rates,E,I,Coords)
rng('shuffle');

%% Phase Attribution Systematic Evaluation
ConnMat = W;

Phases = PhaseDistance(W,Rates);
[~,ix] = sort(Phases);
[~,ixp] = sort(Coords);
Coords = Coords(ixp)';
Rates = Rates(:,ix);
%N.Types = N.Types(ix);
ConnMat = ConnMat(ix,ix);
Orsyn = logical(ConnMat);
orpos = Coords;
%%
Diffs = [];
MeanProjLength =[];
jj = 1;
fig = figure;
for L=1:1:10
    for ii = 1:10
        INT = [];
        for i = 1:L 
            if(i<L)    
                int = orpos((i-1)*round(length(Coords(:))/L)+1:i*round(length(Coords(:))/L),:);
            else
                int = orpos((i-1)*round(length(Coords(:))/L)+1:end,:);
            end
            INT = [INT {int}];
        end
    %    rng(10); % Comment this line out if you want to assign phase aribtrarily. Keep it if you want to reproduce the same phase assignement for each network
        INT = INT(randperm(length(INT)));
        Coords = cell2mat(INT');
        %% Gets the profile of projection for each iteration of the phase assignement
        Diffs(:,ii,jj) = SpatialCoupling(W,E,I,Coords,false);
        %% Computes the mean projection length (metabolic cost) in the network for each iteration of the phase assignement
        DRC = abs(-bsxfun(@minus,Coords(:),Coords(:)'));  
        DRC(~Orsyn) = 0;
        MeanProjLength(:,ii,jj) = median(DRC(DRC>0),'all');
        %ecdf(DRC(DRC>0));
        scatter(L,MeanProjLength(:,ii,jj),'filled');
        hold on
    end
    jj = jj+1;
%box off
%set(gcf, 'WindowState', 'maximized');
%drawnow;  % Ensure the plot is fully rendered
%save2pdf(fig,['./'],['MeanProjLength'],'-dpdf');
end
%%
