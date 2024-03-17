function [SimVector,SystemEigenvector,CovarianceEigenvector,DomTuringMode] = SpatialActivityDistribution(W,Rates)
    
    [V,D] = eig(W);
    EVs = diag(D);
    [~,ind] = sort(real(EVs),'descend');
    Vs = V(:,ind);
    VDom = Vs(:,1);

    %% System eigenvector
    SystemEigenvector = real(VDom);



    %% Principal Component 1
  %  coeffs = pca(Rates);
  %  CovarianceEigenvector = coeffs(:,1);



 %   TuringModes = zeros(101,200);
 %   for k = 0:0.1:10  
    
 %       combsin = zeros(200,1);
 %       combcos = zeros(200,1);
 %       pad = combsin;
 %       combsin(1:1:end) = sin(k*2*pi*(1:length(combsin))/(length(combsin)));%.*(1:-1/(length(combsin)-1):0);
 %       combsin = [pad;combsin;pad];
 %       combcos(1:1:end) = cos(k*2*pi*(1:length(combcos))/(length(combcos)));%.*(1:-1/(length(combcos)-1):0);
 %       combcos = [pad;combcos;pad];
 %       sigsin = movmean(conv(combsin,Diff),50);
 %       sigcos = movmean(conv(combcos,Diff),50);
 %       sigtps = sigsin(((length(Diff)/2)+length(pad)):(length(sigsin)-(length(Diff)/2+length(pad))));
 %       sigtpc = sigcos(((length(Diff)/2)+length(pad)):(length(sigcos)-(length(Diff)/2+length(pad))));

        
 %       k = k*10;
 %       k = cast(k,"int8");
 %       k = k + 1;
 %       TuringModes(k,:) = sigtpc;
    
 %   end
        
    
%    Amps = rms(TuringModes');
%    [~,ix] = max(Amps);
    
    %% Dominant Turing Mode
%    DomTuringMode = TuringModes(ix,:);
%    DomTuringMode = DomTuringMode(51:150);


    %% Spatial distribution from simulation
    SimVector = Rates(5000,:);

    
    fig = figure;
    set(gcf, 'WindowState', 'maximized');
    subplot(2,1,1);
    bar(SystemEigenvector);
%    set(gca,"YDir","reverse");
    title('Eigen Mode','FontSize',16);
    xlabel('1D Spatial Coordinate','FontSize',16);
    ylabel('∝ Activity','FontSize',16);
    box off
    subplot(2,1,2);
    bar(SimVector);
    title('Simulation','FontSize',16);
    xlabel('1D Spatial Coordinate','FontSize',16);
    ylabel('∝ Activity','FontSize',16);
    box off
 
    


%    subplot(2,2,3);
 %   plot(CovarianceEigenvector);
%    set(gca,"YDir","reverse");
 %   title('Dominant Covariance Eigenvector (PC1)','FontSize',16);
 %   xlabel('1D Spatial Coordinate','FontSize',16);
 %   ylabel('∝ Activity','FontSize',16);
 %   box off
 %   subplot(2,2,4);
 %   plot(linspace(0,1000,100),DomTuringMode);
 %   title('Dominant Turing Mode','FontSize',16);
 %   xlabel('1D Spatial Coordinate','FontSize',16);
 %   ylabel('∝ Activity','FontSize',16);
 %   box off

    
%    drawnow;  % Ensure the plot is fully rendered
%    save2pdf(fig,['./'],['SpatialSummary'],'-dpdf');


end

