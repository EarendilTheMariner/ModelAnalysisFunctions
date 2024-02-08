
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


ComplexDoms = [];

%% Save vars in unique iteration folder
for ii=1:100
    clearvars -except ii ComplexDoms
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);
    [W, DominantMode] = BalancedRecRandom('Rates','false');
    save("output.mat","W");
    save("DominantMode.mat","DominantMode");
    if ~isreal(DominantMode)
        ComplexDoms = [ComplexDoms DominantMode];
    else
        disp('not complex');
    end
    cd ..    
    close all
end


%% N=2000 generations


for ii=1:100

    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load("output.mat");
    Rates = R;

    if isfile("Coords.mat")
        Coords = PhaseDistance(W,Rates);
        [SymW,]
        
        end
    else
        delete Coords
        disp("tiny variance, no oscillations")
    end
    cd ..
end


%% Compute kernel and save .png file in each iteration folder
for ii=1:100

    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load("Rates.mat");
    load("output.mat");
    load("MeanVarRates.mat");
    if RatesVar > 1
        Coords = PhaseDistance(W,Rates);        
        if ~isempty(Coords)
            CouplingKernel(W,0.5,Coords,Rates);
            save("Coords","Coords");
            close all
        else
            delete Coords.mat
        end
    else
        disp("No oscillations")
    end
    VarCutOff = 1;
    save("VarCutOff","VarCutOff");
    ProminenceCutoff = 5;
    save("ProminenceCutoff", "ProminenceCutoff");
    cd ..
end


%% Compute number of oscillating networks 
NumOsc_N1000 = 0;
for ii=1:100
    clearvars -except ii NumOsc_N1000 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    cd("BaseNetwork");
    load("MeanRatesVar.mat");
    if MeanRatesVar > 1
        NumOsc_N1000 = NumOsc_N1000 + 1;
    else 
        disp("no oscillations")
    end 
    cd ..
    cd ..
end

%% Prominent peaks and cleanup


for ii = 1:100
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);

    load("Rates.mat")
    [~,scores] = pca(Rates);
    PCRef = scores(:,1);
    RatesVar = mean(var(Rates));
    if RatesVar > 1

        figure;
        findpeaks(PCRef,'MinPeakProminence',5,'Annotate','extents');
        title('Mean unit variance =',num2str(RatesVar));
        set(gcf, 'WindowState', 'maximized');
        drawnow;  % Ensure the plot is fully rendered
        frame = getframe(gcf);  % Get the current frame
        im = frame2im(frame); % Convert the frame to an image
        imwrite(im, 'PeakIdentification.png');
        save("MeanVarRates","RatesVar");
    else
        save("MeanVarRates","RatesVar");
    end
    close all
    cd ..
end

for ii = 1:100
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    if isfile("Kernel.png")
        load("Coords.mat");
        load("output.mat");
        [SymW,SymRates,SymDominantMode] = SymmetricGains(W,Coords);
        save("SymW","SymW");
        save("SymRates","SymRates");  
        save("SymDominantMode","SymDominantMode");
        close all
    else
        disp("no kernel")
    end
    cd ..
end

for ii = 1:100 
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    mkdir("BaseNetwork");
    if isfile("Kernel.png")

        mkdir("SymmetricNetwork");
        mkdir("BalancedSymmetricNetwork");

        movefile output.mat BaseNetwork
        movefile MeanVarRates.mat BaseNetwork
        movefile VarCutOff.mat BaseNetwork
        movefile ProminenceCutoff.mat BaseNetwork
        movefile PeakIdentification.png BaseNetwork
        movefile Kernel.png BaseNetwork
        movefile dynamics.png BaseNetwork
        movefile Coords.mat BaseNetwork
        
        movefile SymRates.mat SymmetricNetwork
        movefile SymW.mat SymmetricNetwork
        movefile SymDominantMode.mat SymmetricNetwork
        movefile SymKernel.png SymmetricNetwork

 %      movefile BalSymRates.mat BalancedSymmetricNetwork
 %      movefile BalSymW.mat BalancedSymmetricNetwork
 %      movefile BalSymDominantMode.mat BalancedSymmetricNetwork
 %      movefile BalSymKernel.png BalancedSymmetricNetwork
    else
        movefile output.mat BaseNetwork
        movefile dynamics.png BaseNetwork
        movefile MeanVarRates.mat BaseNetwork
        movefile VarCutOff.mat BaseNetwork
        movefile ProminenceCutoff.mat BaseNetwork
        if isfile("PeakIdentification.png")
            movefile PeakIdentification.png BaseNetwork
        else
            disp('no peak plot')
        end
    end
    cd ..
end


BaseNetDomReals = 0;
BaseNetDomComp = 0;

SymNetDomReals = 0;
SymNetDomComp = 0;

BalSymNetDomReals = 0;
BalSymNetDomComp = 0;

for ii = 1:100 
    clearvars -except ii BaseNetDomReals BaseNetDomComp SymNetDomReals SymNetDomComp BalSymNetDomReals BalSymNetDomComp
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    if isfolder("SymmetricNetwork\")
        cd BaseNetwork
        load('output.mat')
        eigs = eig(W);
        [~,eix] = max(real(eigs));
        if isreal(eigs(eix))
            BaseNetDomReals = BaseNetDomReals + 1;
        else
            BaseNetDomComp = BaseNetDomComp + 1;
        end
        clearvars -except ii BaseNetDomReals BaseNetDomComp SymNetDomReals SymNetDomComp BalSymNetDomReals BalSymNetDomComp
        cd ..

        cd SymmetricNetwork\
        load("SymDominantMode.mat")
        if isreal(SymDominantMode)
            SymNetDomReals = SymNetDomReals + 1;
        else
            SymNetDomComp = SymNetDomComp + 1;
        end
        clearvars -except ii BaseNetDomReals BaseNetDomComp SymNetDomReals SymNetDomComp BalSymNetDomReals BalSymNetDomComp
        cd ..
        cd BalancedSymmetricNetwork\
        load("BalSymDominantMode.mat")
        if isreal(BalSymDominantMode)
            BalSymNetDomReals = BalSymNetDomReals + 1;
        else 
            BalSymNetDomComp = BalSymNetDomComp + 1;
        end
        cd ..
    end
    cd ..
end
save("ModeCounts");



fig = figure;
BarVals = [1,1,1,1,1,1];
x = categorical({'N = 200', 'N = 500', 'N = 1000', 'N = 2000', 'N = 5000', 'N = 10000'});
x = reordercats(x,string(x));
bar(x,BarVals,0.3);
yLabel = ylabel('f');
title('Fraction Of Complex Dominant Modes')
set(gca, 'xticklabel', x);
set(yLabel, 'FontSize', 18); % Set y-axis label font size
set(gca, 'FontSize', 18);
box off
grid()
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['ComplexFractionDist'],'-dpdf');






ComplexDoms = [];
for ii = 1:100
    clearvars -except ii ComplexDoms
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load("DominantMode.mat");
    if ~isreal(DominantMode)
        ComplexDoms = [ComplexDoms DominantMode];
    else
        disp('Not complex')
    end
    cd .. 
end
save("ComplexDoms.mat","ComplexDoms");



DominantModesN10000 = [];
ComplexCount = 0;
for ii = 1:100
    clearvars -except ii ComplexCount DominantModesN10000
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load("Eigenvalues.mat");
    [~,eix] = max(real(Eigenvalues));
    DominantMode = Eigenvalues(eix);
    save("DominantMode.mat","DominantMode");
    if ~isreal(DominantMode)
        ComplexCount = ComplexCount + 1;
        DominantModesN10000 = [DominantModesN10000 DominantMode];
    else
        disp('real dominant ):')
    end
    
    [brevr,bixev] = sort(real(Eigenvalues),'descend');
    bievr = imag(Eigenvalues(bixev));

    fig = figure('Visible','off');
    scatter(brevr,bievr);
    vline(1);
   % xlim([-1.5,1.5]);
    axis equal
    box off

    set(gcf, 'WindowState', 'maximized');
    drawnow;  % Ensure the plot is fully rendered
    save2pdf(fig,['./'],['Spectrum'],'-dpdf');
    close all
    cd ..
end
save("ComplexCountN10000.mat","ComplexCount");
MeanImagPartN10000 = mean(imag(DominantModesN10000));
save("MeanImagPartN10000.mat","MeanImagPartN10000");



DominantModesN200 = [];
for ii = 1:100
    clearvars -except ii DominantModesN200
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load("DominantMode.mat")
    DominantModesN200 = [DominantModesN200 DominantMode];
    cd .. 
end
save("DominantModesN200.mat","DominantModesN200");
MeanImagPartN200 = mean(imag(DominantModesN200));
save("ImagPartN200.mat","MeanImagPartN200");









figure;
imagesc(W);
currentColormap = colormap;
range = max(W(:)) - min(W(:)); % Range of your data
zeroIndex = -min(W(:)) / range; % Position of zero relative to the data range
zeroColormapIndex = round(zeroIndex * size(currentColormap, 1)) + 1; % Colormap index for zero
currentColormap(zeroColormapIndex, :) = [1, 1, 1]; % Set the color to white
colormap(currentColormap); % Apply the modified colormap


set(gca, 'YDir', 'normal');
xlabel('Neuron j');
ylabel('Neuron i');
title('Network connectivity matrix')
colorbar





GainSweep = [linspace(0.85,1.06,100) linspace(1.06,0.85,100)];
for ii = 1:length(GainSweep)
    close all
    clearvars -except W ii GainSweep 

    gain = ones(size(W,1),20000).*GainSweep(ii);

    Rates = SimulateNetwork(W,20000,'gain',gain);
    Rates = Rates(10001:end,:);

    figure;
    plot(Rates);
    grid();
    xlabel('Time (ms)');
    ylabel('Rate (Hz)');
    xlim([0,10000]);
    ylim([0,50]);
    box off
    if ii <= 100
        legend('Gain ↑')
    else
        legend('Gain ↓')
    end
    vline(1);
    set(gcf, 'WindowState', 'maximized');

    drawnow;  % Ensure the plot is fully rendered
    frame = getframe(gcf);  % Get the current frame
    im = frame2im(frame); % Convert the frame to an image
    imwrite(im, 'temp.png');
    im = imread('temp.png');
    [imind,cm] = rgb2ind(im,256);  % Convert the image to indexed color

    filename = 'GainSweepR.gif';

    % Write to GIF
    if ii == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.1);
    end
end

delete('temp.png');
gif2avi('GainSweepR.gif');


BaseW = W;
GainSweep = [linspace(0.8,1.2,100) linspace(1.2,0.8,100)];
for ii = 1:length(GainSweep)
    close all
    clearvars -except ii GainSweep BaseW

    W = GainSweep(ii).*BaseW;

    bevr = eig(W);
    [brevr,bixev] = sort(real(bevr),'descend');
    bievr = imag(bevr(bixev));

    figure;
    scatter(brevr,bievr);
    vline(1);
    xlim([-1.5,1.5]);
    xlabel('Real');
    ylabel('Imaginary')
    title('Eigenspectrum');
    axis equal tight
    box off
    grid()
   
    if ii <= 100
        legend('Gain ↑')
    else
        legend('Gain ↓')
    end
    vline(1);
    set(gcf, 'WindowState', 'maximized');

    drawnow;  % Ensure the plot is fully rendered
    frame = getframe(gcf);  % Get the current frame
    im = frame2im(frame); % Convert the frame to an image
    imwrite(im, 'temp.png');
    im = imread('temp.png');
    [imind,cm] = rgb2ind(im,256);  % Convert the image to indexed color

    filename = 'GainSweep.gif';

    % Write to GIF
    if ii == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.1);
    end
end

delete('temp.png');
clear





bevr = MeanDomsN;
[brevr,bixev] = sort(real(bevr),'descend');
bievr = imag(bevr(bixev));

figure;
scatter(brevr,bievr);
vline(1);
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5])
xlabel('Real');
ylabel('Imaginary')
title('Eigenspectrum');
axis equal
box off
grid()

DiffMat = [];
for ii = 1:100 
    clearvars -except ii DiffMat
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    cd("BaseNetwork");
    if isfile("Diff.mat")
        load('Diff.mat');
        DiffMat = [DiffMat, difference];

    else
        disp("No oscillations")
    end
    cd .. 
    cd ..


end
dataPos = zeros(1,200);
PosIx = MeanDiff > 0;
dataPos(PosIx) = MeanDiff(PosIx);
fig = figure(Position=screensize);
bar([-1000:10:990],dataPos,'FaceColor',[Colors().BergElectricBlue],'FaceAlpha',0.5);
hold on

dataNeg= zeros(1,200);
NegIx = MeanDiff < 0;
dataNeg(NegIx) = MeanDiff(NegIx);
bar([-1000:10:990], dataNeg,'FaceColor',[Colors().BergOrange],'FaceAlpha',0.5);
legend({'Excitation','Inhibition'});




screensize = get(groot,'ScreenSize');
edges = [-1000:10:1000];


clearvars -except MeanDiff
for ii = 1:100
    clearvars -except ii MeanDiff
    folderName = sprintf('Iteration_%d', ii);

    cd(folderName);
    [W, Diff, Eigenvalues, DominantMode, Coords, I, E] = BalancedRecStructured(500,MeanDiff);
    save('W_2.mat',"W");
    save("Diff_2.mat","Diff");
    save("Eigenvalues_2.mat","Eigenvalues");
    save("DominantMode_2.mat","DominantMode");
    save("Coords_2.mat","Coords");
    save("I_2.mat","I");
    save("E_2.mat","E");
    cd ..
end


figure;
Eigenvalues = eig(W);
[brevr,bixev] = sort(real(Eigenvalues),'descend');
bievr = imag(Eigenvalues(bixev));

figure;
scatter(brevr,bievr);
vline(1);
% xlim([-1.5,1.5]);
axis equal
box off


GainSweep = linspace(0,3,100);
for ii = 1:length(GainSweep)
    close all
    clearvars -except ii Coords E I W GainSweep
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);    


    [WCtrl,DiffCtrl,EigenvaluesCtrl,RatesCtrl] = BumpControl(W,Coords,I,E,1,GainSweep(ii));
    save("WCtrl.mat", "WCtrl")
    save("DiffCtrl.mat","DiffCtrl");
    save("RatesCtrl.mat","RatesCtrl");
    save("EigenvaluesCtrl.mat","EigenvaluesCtrl");
    BumpScaling = GainSweep(ii);
    save("BumpScaling.mat","BumpScaling")
    cd ..
end
    

for ii = 1:99
    close all
    clearvars -except ii
    folderName = sprintf('Iteration_%d', ii+1);
    cd(folderName);    
    load("DiffCtrl.mat");
    load("EigenvaluesCtrl.mat");
    load("RatesCtrl.mat");

    cd ..

    [brevr,bixev] = sort(real(EigenvaluesCtrl),'descend');
    bievr = imag(EigenvaluesCtrl(bixev));

    % Visual overview
    
    fig = figure();
    
    subplot(2,2,1);
    x = linspace(-1000,990,200);
    bar(x,DiffCtrl);
    ylim([-30 70])
    box off
    
    subplot(2,2,2);
    scatter(brevr,bievr);
    vline(1);
    % xlim([-1.5,1.5]);
    axis equal
    box off
    
    subplot(2,2,3:4);
    plot(RatesCtrl);
    grid();
    xlabel('Time (ms)');
    ylabel('Rate (Hz)');
    xlim([0,10000]);
    ylim([0,70]);
    box off

    set(gcf, 'WindowState', 'maximized');
    drawnow;  % Ensure the plot is fully rendered
    frame = getframe(gcf);  % Get the current frame
    im = frame2im(frame); % Convert the frame to an image
    imwrite(im, 'temp.png');
    im = imread('temp.png');
    [imind,cm] = rgb2ind(im,256);  % Convert the image to indexed color

    filename = 'Bump4Sweep.gif';

    % Write to GIF
    if ii == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.2);
    end
end




I = logical(I);
E = logical(E);
WRaw = W;
Dist = bsxfun(@minus,Coords,Coords');
L = range(Coords,2);
edges = -L:10:L;
InhNegBump = (Dist(:,:) < 0) & (WRaw(:,:) < 0);


BumpScale = linspace(0,3,100);
for ii = 1:100
    close all
    clearvars -except ii WRaw E I BumpScale InhNegBump
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);   

    W = WRaw;

    W(InhNegBump) = BumpScale(ii).*W(InhNegBump);

    W = BalanceConnectivity(W);

    Eigenvalues = eig(W);

    [~,eix] = max(real(Eigenvalues));
    
    DominantMode = Eigenvalues(eix);
    
    GlobalScaling = 1/(real(DominantMode)*0.95);
    save("GlobalScaling.mat","GlobalScaling");
    BumpScaling = BumpScale(ii);
    save("BumpScaling.mat","BumpScaling");
    cd ..
end

BumpScaleFactors = [];
GlobalScaleFactors = [];
for ii = 1:100
    clearvars -except ii BumpScaleFactors GlobalScaleFactors
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);   
    load("BumpScaling.mat");
    load("GlobalScaling.mat");
    BumpScaleFactors = [BumpScaleFactors BumpScaling];
    GlobalScaleFactors = [GlobalScaleFactors GlobalScaling];
    cd ..
end




figure;
scatter(BumpScaleFactors,GlobalScaleFactors);
ylabel("Global Scaling Factor (after balancing)")
xlabel("Bump Scaling Factor")



