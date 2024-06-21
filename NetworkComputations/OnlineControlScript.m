

t = 20000;
Input = ones(size(W,1),t).*20;
GlobalSweep = 10:1:40;
for ii=1:length(GlobalSweep)
    close all
    clearvars -except ii GlobalSweep W PopMat t Input
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);






 %   Input(logical(E),:) = 10;
%    Input(logical(I),:) = 0;

%    Input(logical(I),10000:22000) = 20;
%    Input(logical(I),22000:35000) = 40;


%
    Input(logical(PopMat(:,2)),1:10000) = 0;
    Input(logical(PopMat(:,2)),10001:end) = 30;
    Input(logical(PopMat(:,2)),20001:30000) = 15;
    Input(logical(PopMat(:,2)),30001:40000) = 30;


    

%    Input(logical(PopMat(:,:)),30000:end) = 40;


%    Input(:,1:10000) = 1;
%    Input(:,10000:20000-1) = 1.2;

%    Input(:,20000:30000-1) = 1.4;

%    Input(:,30000:40000-1) = 1.6;
%    Input = repmat(Input,1000,1);


    


   % Input(logical(PopMat(:,3)),15000:end) = Sp3Sweep(ii);
    
    I_e = GlobalSweep(ii).*Input;
    Rates = SimulateNetwork(W.*1.5,t,'I_e',Input);

    
    save("Rates.mat","Rates");
    save("Input.mat","I_e");
    
    fig = figure;
    subplot(2,1,1);
    plot(Rates);
    ylim([0 70]);
    box off
    
    subplot(2,1,2);
    plot(Input(2,:));
    ylim([0 70]);
    hold on
    plot(Input(1,:));
    %box off
    %plot(InputDown(3,:));
    %legend('Sp2 Input', 'Sp4 Input','Sp3 Input');
    legend('Global Input');
 
    fig = figure;
    set(gcf, 'WindowState','maximized');
    box off
    plot(Input(3,10001:30000));
    ylim([0 40]);
    drawnow;  % Ensure the plot is fully rendered
    save2pdf(fig,['./'],['Input'],'-dpdf');
    
    cd ..
end


clear
load('BaseNetwork.mat');
t = 20000;
Input = ones(size(W,1),t).*20;
Sp3Sweep = 0:1:20;
for ii=1:length(Sp3Sweep)
    close all
    clearvars -except ii Sp3Sweep W Input PopMat t
    folderName = sprintf('Iteration_%d', ii);
    mkdir(folderName);
    cd(folderName);


    %Input(logical(E),:) = 60;
    Input(logical(PopMat(:,3)),1:end) = Sp3Sweep(ii);
    %Input(logical(PopMat(:,1)),1:end) = Sp3Sweep(ii);
   % Input(logical(PopMat(:,5)),1:end) = Sp3Sweep(ii);
    %Input(logical(PopMat(:,4)),15000:end) = 25;
    %Input(logical(PopMat(:,3)),15000:end) = 40;
    
    
    Rates = SimulateNetwork(W.*1.5,t,'I_e',Input);
    save("Rates.mat","Rates");
    save("Input.mat","Input");
    
    
    fig = figure;
  %  subplot(2,1,1);
    plot(Rates)
    ylim([0 70]);
    box off
    
 %   subplot(2,1,2);
 %   plot(Input(3,:));
 %   ylim([-10 50]);
 %   hold on
 %   plot(Input(1,:));
 %   box off
    %plot(InputDown(3,:));
    %legend('Sp2 Input', 'Sp4 Input','Sp3 Input');
 %   legend('Sp4 Input','Global Input');
    
    
    set(gcf, 'WindowState','maximized');
    drawnow;  % Ensure the plot is fully rendered
    save2pdf(fig,['./'],['OnlineFreqSweep'],'-dpdf');
    
    cd ..
end

for ii=1:91
    close all
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load('Rates');

    [~,scores] = pca(Rates(1000:14999,:));
    FrequencyStart = FrequencyFromFFT(scores(:,1));
    save('FrequencyStart.mat','FrequencyStart');

    [~,scores] = pca(Rates(15000:end,:));

    if var(scores(15000:end,1)) > 1
        FrequencyEnd = FrequencyFromFFT(scores(:,1));
        save('FrequencyEnd.mat','FrequencyEnd');
    else 
        FrequencyEnd = NaN;
        save('FrequencyEnd.mat','FrequencyEnd');
    end
    cd ..
end



Freqs = zeros(1,41);

for ii=1:41
    close all
    clearvars -except Freqs ii
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load('Frequency');
    Freqs(ii) = Frequency;
    cd ..
end

for ii=1:1:41
    %close all
    clearvars -except ii 
    folderName = sprintf('Iteration_%d', ii);
    cd(folderName);
    load('Rates');

    RatesStart = Rates(1000:end,:);
    [~, Frequency] = ComputeFrequency(RatesStart);
    save('Frequency.mat',"Frequency");
    cd ..
end


fig = figure; 
scatter(0:40,AsymFreqs/max(AsymFreqs),'filled'); ylim([0 1]); xlim([0 40]); 
%xlabel('Sp2 Step Input Change'); ylabel('Relative Frequency Change');
hold on 
scatter(0:40,SymFreqs/max(SymFreqs),'filled'); ylim([0 1]); xlim([0 40]); 
%xlabel('Sp2 Step Input Change'); ylabel('Relative Frequency Change');
%legend('Population Asymmetry', 'Population Symmetry');
set(gcf, 'WindowState','maximized');
drawnow;  % Ensure the plot is fully rendered
save2pdf(fig,['./'],['FreqVsSp2Input'],'-dpdf');

a_1 = StepInputChange;
a_2 = RC;

b_1 = StepInputChange;

b_2 = RC;







