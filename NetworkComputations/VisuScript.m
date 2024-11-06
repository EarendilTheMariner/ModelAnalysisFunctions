
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
    clearvars -except ii Mean
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
BarVals = [Bal, NoBal];
x = categorical({'Balanced', 'Unbalanced'});
x = reordercats(x,string(x));
bar(x,BarVals,0.1);
yLabel = ylabel('Norm. Oscillatory Gain Range');
yticks([0 0.5 1]);
%title('Fraction Of Complex Dominant Modes')
set(gca, 'xticklabel', x);
set(yLabel, 'FontSize', 18); % Set y-axis label font size
set(gca, 'FontSize', 18);
colororder("earth")
ylim([0 1]);
box off
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['OscillatoryRange'],'-dpdf');






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

%% BUMP CONTROL LOOP
PopMat = logical(PopMat);
GainSweep = linspace(0.1,2,100);
for ii = 1:length(GainSweep)
    close all
    clearvars -except ii Coords E I W GainSweep PopMat 
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);
    WMod = W
    WMod = W(:,PopMat(:,4)).*GainSweep(ii);



    cd ..
end



cd ..
cd BumpSweep2&4NoBal

GainSweep = linspace(0.1,2,100);
for ii = 1:length(GainSweep)
    close all
    clearvars -except ii Coords E I W GainSweep 
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);    


    [WCtrl,DiffCtrl,EigenvaluesCtrl,RatesCtrl] = BumpControl(W,Coords,I,E,2,GainSweep(ii));
    save("WCtrl.mat", "WCtrl")
    save("DiffCtrl.mat","DiffCtrl");
    save("RatesCtrl.mat","RatesCtrl");
    save("EigenvaluesCtrl.mat","EigenvaluesCtrl");
    BumpScaling = GainSweep(ii);
    save("BumpScaling.mat","BumpScaling")
    cd ..
end



%% GIF LOOP
for ii = 1:185
    close all
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);    
 %   load("DiffCtrl.mat");
 %   load("EigenvaluesCtrl.mat");
    load('Output.mat');
  %  load("Rates.mat");

    cd ..

    % Visual overview
    
    fig = figure;
    
  %  ha = axes('units','normalized','position',[0 0 1 1]);
  %  uistack(ha,'bottom');
   % imshow('Background_V2.png');
   % truesize;
  %  set(ha,'handlevisibility','off','visible','off');
    fig = gcf;
    
    subplot(1,3,2);
    SpatialCoupling(W,E,I,Coords,true);
    set(gca, 'Color', 'none');
    ylim([-20 50]);
    yticks([-30 0 30 60]);
    box off;
 
    
    subplot(1,3,3);
    EigenSpectrum(W);
    set(gca, 'Color', 'none');
    xticks([0 0.5 1]);
    yticks([-0.25 0 0.25]);
    box off;
    
    subplot(1,3,1);
    plot(PDFs);
    box off;
    
    set(fig, 'WindowState', 'maximized');
    
    % Save the figure directly to a file
    print(fig, 'temp.png', '-dpng', '-r300');  % Save at 300 dpi
    
    im = imread('temp.png');
    [imind, cm] = rgb2ind(im, 256);  % Convert the image to indexed color
    
    filename = 'Overview.gif';
    
    % Write to GIF
    if ii == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end

%% FREQUENCY LOOP
Amps = [];
Freqs = [];
BumpScalars = [];

for ii = 1:100
  %  close all
    clearvars -except ii Freqs BumpScalars 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);

    load("Rates.mat");
    load("BS2.mat");
    

    [~,Frequency] = ComputeFrequency(Rates(10000:end,:));


    Freqs = [Freqs Frequency];
    BumpScalars = [BumpScalars BS2];

    cd ..
end




%% FREQUENCY SUMMARYS FIGURE
fig = figure;
scatter(linspace(0.01,3,100),NormImagParts/(2*pi));
hold on
scatter(linspace(0.01,3,100),Freqs,[]);
legend('Imaginary Part','∼ Frequency','Location','northeast');
legend('boxoff');
legend('FontSize',20);
%ylabel('Norm. Units','FontSize',20);
xlabel('Bump 4 Scale Factor', 'FontSize',20);
%ylabel('Norm. Units','FontSize',20);
box off

set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['RawFrequencySummary'],'-dpdf');

%% SPATIAL ACTIVITY LOOP

Amplitudes = zeros(41,1);
for ii = 1:41
    clearvars -except ii Amplitudes
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load('Rates.mat');
    Rates = Rates(5001:end,:);
    MR = mean(Rates);
    MR = repmat(MR,10000,1);
    Rates = Rates-MR;
    Amplitude = rms(Rates);
    Amplitude = mean(Amplitude);

    Amplitudes(ii) = Amplitude;

  %  [~, Frequency] = ComputeFrequency(RatesCtrl);
  %  Freqs(ii,1) = Frequency;
    cd ..
end

%% SP2 vs Sizes
fig = figure;
set(gcf, 'WindowState', 'maximized');
plot(linspace(0.1,2,100),N200Amps/max(N200Amps),'k');
hold on
plot(linspace(0.1,2,100),N1000Amps/max(N1000Amps),'k');
hold on 
plot(linspace(0.1,2,100),N200Freqs/max(N200Freqs),'r');
hold on
plot(linspace(0.1,2,100),N1000Freqs/max(N1000Freqs),'r');
hold on 
plot(linspace(0.5,2,30),N5000Amps/max(N5000Amps),'k');
hold on
plot(linspace(0.5,2,30),N5000Freqs/max(N5000Freqs),'r');
%legend('N = 200 Amplitude','N = 1000 Amplitude','N = 200 Frequency','N = 1000 Frequency','N = 5000 Amplitude','N = 5000 Frequency','Location','northwest');
xlabel('Sp2 Gain');
ylabel('Norm. Frequency/Amplitude');
xlim([0 2]);
ylim([0 1]);
box off
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['N200_N1000__N5000_Freq_Amplitude_Sp3_Curves'],'-dpdf');

%% SP3 vs Sizes
fig = figure;
set(gcf, 'WindowState', 'maximized');
scatter(linspace(0.1,2,100),N200Amps/max(N200Amps),[],'black','filled');
hold on
scatter(linspace(0.1,2,100),N1000Amps/max(N1000Amps),"^",'black','filled');
hold on 
scatter(linspace(0.1,2,100),N200Freqs/max(N200Freqs),[],'red','filled');
hold on
scatter(linspace(0.1,2,100),N1000Freqs/max(N1000Freqs),"^",'red','filled');
hold on 
scatter(linspace(0.5,2,30),N5000Amps/max(N5000Amps),'square','black','filled');
hold on
scatter(linspace(0.5,2,30),N5000Freqs/max(N5000Freqs),'square','red','filled');
legend('N = 200 Amplitude','N = 1000 Amplitude','N = 200 Frequency','N = 1000 Frequency','N = 5000 Amplitude','N = 5000 Frequency','Location','northwest');
xlabel('S3 Gain');
ylabel('Norm. Frequency/Amplitude');
xlim([0 2]);
ylim([0 1]);

drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['NN200_N1000__N5000_Freq_Amplitude_Sp3_'],'-dpdf');



%% SPATIAL DISTRIBUTION VS FREQUENCY PLOTS
fig = figure;
set(gcf, 'WindowState', 'maximized');
Sweep = round(linspace(7,100,9));
for ii = 1:length(Sweep)
    clearvars -except ii fig Sweep ImagParts
    folderName = sprintf('Iteration_%d', Sweep(ii));
    cd(folderName);

    load('RatesCtrl.mat');
    load('BumpScaling.mat');

    if ii == 1

        subplot(3,3,ii)
        imagesc(flipud(RatesCtrl));
        box off
        set(gca, 'YDir', 'normal'); 
        title(['Bump Scaling Factor = ', num2str(round(BumpScaling,2)), ', Norm. Frequency = ', num2str(round(ImagParts(Sweep(ii)),2))],'FontSize',16);
        xlabel('1D Spatial Coordinate','FontSize',16);
        ylabel('Time (ms)','FontSize',16);
        yticks([1 5000 10000]);
        box off
        hold on
        cd ..  
    else 
        subplot(3,3,ii)
        imagesc(flipud(RatesCtrl));
        box off
        title(['Bump Scaling Factor = ', num2str(round(BumpScaling,2)), ', Norm. Frequency = ', num2str(round(ImagParts(Sweep(ii)),2))],'FontSize',16);
        yticks([1 5000 10000]);
        set(gca, 'YDir', 'normal');
        hold on
        cd ..
    end
end
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['SpatialActivitySampling'],'-dpdf');








fig = figure;

scatter(linspace(0.001,2,100),AmpsCenter);
xlim([0 2]);
ylim([20 45]);
ylabel('Oscillation Amplitude');
xlabel('sp3 Gain')
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['AmplitudeVsSP3'],'-dpdf');



fig = figure;
scatter(linspace(0.001,2,100),FreqsCenter);
xlim([0 2]);
ylim([0 1]);
ylabel('Oscillation Frequency');
xlabel('sp3 Gain');
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['Frequency&SP3'],'-dpdf');

fig = figure;
scatter(linspace(0.001,2,100),AmpsB4);
xlim([0 2]);
ylim([20 45]);
ylabel('Oscillation Amplitude');
xlabel('sp4 Gain');
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['AmplitudeVsSP4'],'-dpdf');



fig = figure;
scatter(linspace(0.001,2,100),Freqs);
xlim([0 2]);
ylim([0 1]);
ylabel('Oscillation Frequency');
xlabel('sp4 Gain');


set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['FrequencyVsSP4'],'-dpdf');


BS = [];
for ii = 1:100
    clearvars -except ii BS
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);

    load('RatesCtrl.mat');
    [~,Frequency] = ComputeFrequency(RatesCtrl);
    if ~isnan(Frequency)
        load('BumpScaling.mat');
        BS = [BS BumpScaling];
    else 
        disp('no oscillations')
    end
   
    cd ..
end


%% Downsampled plot
downsampled = downsample(RatesCtrl,50);
fig = figure; 
plot(linspace(1,10000,200),downsampled);
xlabel('Time','FontSize',20);
ylabel('Rate','FontSize',20);
xlim([0,10000]);
ylim([0,70]);
box off
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['DownsampledRates_32'],'-dpdf');




%% Projectome vs. Network Length
ExtensionSweep = 1:1:10;
for ii = 1:length(ExtensionSweep)
    clearvars -except ii MeanDiff ExtensionSweep
    close all

    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);

    sp4 = MeanDiff(123:185);

    sp4_ext = kron(sp4,ones(ExtensionSweep(ii),1));

    PDF = [zeros(1,length(sp4_ext)-length(sp4)) MeanDiff(1:122)' sp4_ext' MeanDiff(186:end)'];

    [W, Diff, Eigenvalues, Coords, I, E] = BalancedRecStructured(8000,PDF',8000);

    EigenStructures = EigenStructure(W,5);
    save("EigenStructures.mat","EigenStructures");
    save("W.mat","W");
    save("Diff.mat","Diff");
  %  save("Rates.mat","Rates");
    save("Coords.mat","Coords");
    save("I.mat","I");
    save("E.mat","E");
    save("Eigenvalues.mat","Eigenvalues");
    EigenSpectrum(W);

    fig = figure; 
    subplot(2,1,1);
    SpatialCoupling(W,E,I,Coords);
    subplot(2,1,2);
    bar(EigenStructures(:,1));
    set(gcf, 'WindowState', 'maximized');
    drawnow;  % Ensure the plot is fully rendered
    save2pdf(fig,['./'],['SpatialOverview'],'-dpdf');



    %SpatialActivityDistribution(W,Rates,Diff);
    cd ..
end
    
cd ..
cd ProjectomeSweep/


ProjectomeSweep = [2/10 4/10 6/10 8/10 10/10 12/10 14/10 16/10 18/10 20/10];

for ii = 1:length(ProjectomeSweep)
    clearvars -except ii MeanDiff ProjectomeSweep LengthSweep
    close all
    
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);
    PDF = kron(MeanDiff,ones(ProjectomeSweep(ii),1));
    [W, Diff, Eigenvalues, Rates, Coords, I, E] = BalancedRecStructured(10000,PDF,10000);
    save("W.mat","W");
    save("Diff.mat","Diff");
    save("Rates.mat","Rates");
    save("Coords.mat","Coords");
    save("I.mat","I");
    save("E.mat","E");
    save("Eigenvalues.mat","Eigenvalues");

    SpatialActivityDistribution(W,Rates,Diff);
    cd ..
end



for ii = 1:10
    clearvars -except ii
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load('RatesCtrl.mat');
    coeffs = pca(RatesCtrl);
    PC1 = coeffs(:,1);
    save('PC1.mat','PC1');
    cd ..
end





SpatialFrequency = zeros(1,10);
for ii = 1:10
    clearvars -except ii SpatialFrequency
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load('PC1.mat');
    Frequency = FrequencyFromFFT(PC1);
    SpatialFrequency(ii) = Frequency;
    cd ..
end

Ratios = [2/10 4/10 6/10 8/10 10/10 12/10 14/10 16/10 18/10 20/10];

for ii = 1:10
    clearvars -except ii Ratios
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load('W.mat');
    
    cd ..
end






fig = figure;
subplot(3,1,1);
bar(downsample_v2);
ylim([-0.035 0.035]);
box off
subplot(3,1,2);
bar(downsample_v5);
ylim([-0.035 0.035]);
box off
subplot(3,1,3);
bar(downsample_v8);
ylim([-0.035 0.035]);
box off 
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['EigenStructer_Iter_2_5_8'],'-dpdf');






load('RatesCtrl.mat');
fig = figure;
plot(RatesCtrl(:,500));
xlim([0 10000]);
ylim([0 70]);
box off 
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['Neuron500Example'],'-dpdf');


fig = figure;
scatter([LowerBumpScalingNoBal/2, UpperBumpScalingNoBal/2], [1 1],'filled');
xlim([0 1]);
hold on 
scatter([LowerBumpScalingBal/2, UpperBumpScalingBal/2], [2 2],'filled');
xlim([0 1]);
box off
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['OscillatoryRanges'],'-dpdf');


%% Paper Figure 3 Sweeps

% Size vs Mean Imag Stats
for ii = 1:100
    close all
    clearvars -except ii MeanDiff
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);
    [W,Coords,Diff,PopMat,E,I] = PopulationSplitNetwork(10000,1000,MeanDiff);
    save('Output.mat');
    DominantMode = eigs(W,1);
    save('DominantMode.mat','DominantMode');
    cd ..
end

ComplexCount = 0;
for ii = 1:100
    clearvars -except ii ComplexCount
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load('DominantMode.mat');
    if ~isreal(DominantMode)
        ComplexCount = ComplexCount + 1;
    else
        disp('non-oscillatory')
    end
    cd ..
end

save("ComplexCount.mat","ComplexCount");


ImaginaryParts = zeros(1,100);
for ii = 1:100
    clearvars -except ii ImaginaryParts
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load('DominantMode.mat');

    ImaginaryParts(ii) = imag(DominantMode);
    cd ..
end
save("ImaginaryParts.mat","ImaginaryParts");
MeanImaginaryParts = mean(ImaginaryParts);
save("MeanImaginryParts.mat","MeanImaginaryParts");


    
fig = figure;
bar(1:6,Frequencies/max(Frequencies));
box off
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['NormExpectedFrequency'],'-dpdf');

% Offline Control


S2 = abs(sum(W(:,logical(PopMat(:,2))),'all'));
S4 = abs(sum(W(:,logical(PopMat(:,4))),'all'));
STot = S2+S4;
Sp4Sweep = linspace(0.1,1.85,100);
for ii = 1:100
    close all
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
 %   mkdir(folderName);
    cd(folderName);   
    load('Output.mat');
%   WMod = W.*1.1;

%   BS4 = Sp4Sweep(ii);
 %  WMod(:,logical(PopMat(:,4))) = WMod(:,logical(PopMat(:,4))).*BS4;
 
 %  BS2 = (STot-BS4*S4)/S2;
 %  WMod(:,logical(PopMat(:,2))) = WMod(:,logical(PopMat(:,2))).*BS2;
%   WMod = BalanceConnectivity(WMod);

    Rates = SimulateNetwork(WMod,20000);
    save("Output.mat","WMod","Rates");

    cd ..
end


Sp3Sweep = linspace(0.5,2,30);
for ii = 1:100
    close all
    clearvars -except ii WRef PopMat E I Coords Sp3Sweep
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);   

    W = WRef;


    W(:,logical(PopMat(:,3))) = W(:,logical(PopMat(:,3))).*Sp3Sweep(ii);

    BS = Sp3Sweep(ii);
    save("BS.mat","BS");

    W = BalanceConnectivity(W);

    save("W.mat","W");

    Rates = SimulateNetwork(W,20000);
    save("Rates.mat","Rates");

    cd ..
end

SizeSweep = 1000:1000:10000;
for ii = 1:10
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
  %  mkdir(folderName);
    cd(folderName);
    load("Output.mat");
    Rates = SimulateNetwork(W,15000);
    save("Rates.mat","Rates");
    

 %   [W,Coords,Diff,PopMat,E,I] = PopulationSplitNetwork(SizeSweep(ii),SizeSweep(ii),MeanDiff);
 %   save("Output.mat");    
     
    cd ..
end

WaveCount = zeros(1,10);
for ii = 1:10
    clearvars -except ii WaveCount
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);
    load("Rates.mat");
    Frequency = FrequencyFromFFT(Rates(end,:));
    WaveCount(ii) = Frequency*length(Rates(end,:));
    cd ..
end


%% Downsampled plot
downsampled = downsample(downsampled',10);
fig = figure;
imagesc(downsampled);
colormap(gray);
%bar(downsampled);
%ylim([0 70]);
box off
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['ActivityVsSpace'],'-dpdf');

    


downsampled = downsample(V,4);
fig = figure;
%imagesc(Rates');
%colormap(gray)
bar(downsampled);
ylim([-0.035 0.035]);
box off


fig = figure;
set(gcf, 'WindowState', 'maximized');
[~,scores] = pca(Rates(10001:end,:));
plot(scores(:,1),scores(:,2));
box off
axis equal
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['PCAExample'],'-dpdf');

    


fig = figure;
bar(N.Projectome);
vline(101)



box off
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['Iter_50_SkewIllustration'],'-dpdf');





comp_count = 0;
for ii=1:100

% Define the range of x-values
x = linspace(-5, 5, 100);

% Compute the sigmoid function
sigmoid = 1 ./ (1 + exp(-x));

% Compute the tanh function
tanh_values = tanh(x);

% Create the plot
fig = figure; % Opens a new figure window
hold on; % Holds the plot for multiple plots on the same graph
plot(x, sigmoid, 'b-', 'LineWidth', 2); % Plot sigmoid in blue with a line width of 2
plot(x, tanh_values, 'r--', 'LineWidth', 2); % Plot tanh in red dashed line with a line width of 2
hold off;
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
% Adding labels and title
xlabel('x');
ylabel('Function value');
title('Comparison of tanh() and Sigmoid Functions');
legend('Sigmoid', 'tanh', 'Location', 'best');
save2pdf(fig,['./'],['ActivationFunctions'],'-dpdf');


% Parameters
fig = figure;
sigma_exc = 1.5;      % Standard deviation of the excitation Gaussian
sigma_inh = 2;        % Standard deviation of the inhibition Gaussian
mu_exc = 0;           % Mean of the excitation Gaussian (biased)
mu_inh = 0;           % Mean of the inhibition Gaussian
a_exc = 2;            % Amplitude of excitation
a_inh = -1;           % Amplitude of inhibition

% Define the range of x values
x = -5:0.1:5;

% Calculate the Gaussian functions
gaussian_exc = a_exc * exp(-(x - mu_exc).^2 / (2 * sigma_exc^2));
gaussian_inh = a_inh * exp(-(x - mu_inh).^2 / (2 * sigma_inh^2));

% Define the symmetric mexican hat function as the sum of excitation and inhibition
mexican_hat = gaussian_exc + gaussian_inh;

% Plot the result
plot(x, mexican_hat, 'LineWidth', 2);
xlabel('x');
ylabel('Amplitude');
title('Symmetric Mexican Hat Function with Excitation and Inhibition');
grid on;
box off
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['AsymMexicanHat'],'-dpdf');


fig = figure;
plot(BumpScalars, BalFreqs/max(BalFreqs));
hold on
plot(BumpScalars, NoBalFreqs/max(BalFreqs));
hold on
xlim([0 2]);
ylim ([0 1]);
box off
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['BalanceConditionDynamics'],'-dpdf');


fig = figure;
bar(N.Projectome);
ylim([-40 50]);
box off
set(gcf, 'WindowState', 'maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['Projectome_79'],'-dpdf');


fig = figure;
plot(ints,Skews/max(abs(Skews))); hold on; plot(ints,[-Freqs(1:50) Freqs(51:end)]/max(abs(Freqs)))
set(gcf, 'WindowState', 'maximized');
box off
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['Freq&Skew_Vs_Sp3Shift'],'-dpdf');




% Define the custom colormap
% Create a colormap with white at zero, black for negative, and dark rose for positive values
negColors = [linspace(0, 1, 128)', linspace(0, 1, 128)', linspace(0, 1, 128)']; % From black to white
darkRoseColor = [0.8, 0, 0.4]; % RGB for dark rose
posColors = [linspace(1, darkRoseColor(1), 128)', linspace(1, darkRoseColor(2), 128)', linspace(1, darkRoseColor(3), 128)']; % From white to dark rose
customColormap = [negColors; posColors];

% Apply the custom colormap
colormap(customColormap);

% Adjust the color limits to be symmetric around zero
maxAbsValue = max(abs(W(:)));
caxis([-maxAbsValue maxAbsValue]);

% Add colorbar for reference
colorbar;


for jj=1:15
    folderName = sprintf('Iteration_%d', jj);  
    cd(folderName);
    load("Projectome.mat");
    ProjectomeOne = Projectome(:,1);
    Sweep = 1:20;
    for ii=1:length(Sweep)
        clearvars -except ii Sweep ProjectomeOne
        folderName = sprintf('Iteration_%d', ii);
        mkdir(folderName);
        cd(folderName);
    %    Projectome = Diffs(:,:,ii);
    %    save("Projectome.mat","Projectome");
    %    MeanAxonLength = MeanProjLength(:,:,ii)
    %    save("MeanAxonLength.mat","MeanAxonLength");
    %    load("Projectome.mat","Projectome");
        [W,Coords,Diff,E,I] = BalancedRecStructured(1000,1000,ProjectomeOne);
        save("Ouput");
        %Rates = SimulateNetwork(W,20000);
        %save("Rates.mat","Rates");
        cd ..
    end
    cd ..
end

Modes = [];
SymIx = [];

Freqs = [];
Sweep = 1:185;
for ii=1:length(Sweep)
    clearvars -except ii Sweep Freqs
    folderName = sprintf('Iteration_%d', ii);
   % mkdir(folderName);
    cd(folderName);
  %  [W,Coords,Diff,PopMat,E,I,PDFs] = SymmetricNetwork(1000,1000,MeanDiff,Sweep(ii));
   % save("Output.mat");
    %DominantMode = eigs(W,1);
 %   save("DominantMode.mat","DominantMode");
 %   load("Output.mat");
 %   R = SimulateNetwork(W,18000);
 %   save("Rates.mat","R");
    load("DominantMode.mat");
    Frequency = abs(imag(DominantMode));
%    load("Rates.mat");
%    R = R(8001:end,:);
%    [Per, Frequency] = ComputeFrequency(R)

    Freqs = [Freqs Frequency];
   % save("Output.mat");
  %  R = SimulateNetwork(W,20000);
  %  save("Rates.mat","R");

    cd ..
end

Skews = [];
for ii=1:100
    clearvars -except ii Coords E I Skews
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load("W.mat");
    Skew = SkewnessFactor(W,E,I,Coords);
    Skews = [Skews Skew];

    cd ..
end

fig = figure; 
plot([1:2:40],PermMeans);

ylim([160 300]);
xlim([0 40]);
set(gcf, 'WindowState', 'maximized');
box off
drawnow
save2pdf(fig,['./'],['MeanOfPermutations'],'-dpdf');




