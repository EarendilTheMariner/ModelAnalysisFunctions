
function [W, DominantMode] = BalancedRecRandom(varargin)
    close all
    rng('shuffle');
    gain = 1;
    Var = 0.1;
    Prob = 0.1;
    n = 200;
    f_exc = 0.5;
    Rates = true;

    % Parsing variable arguments 
    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'Gain'
                gain = varargin{ii+1};
            case 'Var'
                Var = varargin{ii+1};
            case 'Prob'
                Prob = varargin{ii+1};
            case 'Size'
                n = varargin{ii+1};
            case 'F_exc'
                f_exc = varargin{ii+1};
            case 'Rates'
                Rates = varargin{ii+1};
        end
    end
    
    % Transmitter arrays


    fE = f_exc*n;
    fI = (1-f_exc)*n;
    E = false(1,n);
    I = E;
    E(1:fE) = true;
    I(fE+1:end) = true;

    % Connectivity matrix and sparsity
    ConProb = ones(n,n)*Prob;
    Sparsifier = binornd(1,ConProb);
    Connectivity = normrnd(gain,Var,size(ConProb)); 
    Connectivity(:,E) = abs(Connectivity(:,E));
    Connectivity(:,I) = -(Connectivity(:,I)*1.1);
    
    Connectivity = Connectivity.*Sparsifier;
    Connectivity = Connectivity - diag(diag(Connectivity)); 
    Connectivity = BalanceConnectivity(Connectivity);
    Sparsity = nnz(~Connectivity)/numel(Connectivity);
    W = Connectivity;
    
    if Rates == true
        % Simulation
        Rates = SimulateNetwork(W,20000);
        Rates = Rates(10001:end,:);
    
    
        % Eigenspectrum
        bevr = eig(W);
        [brevr,bixev] = sort(real(bevr),'descend');
        bievr = imag(bevr(bixev));
    
        % PCA
        [~, scores] = pca(Rates);
       
        % Plots 
        figure;
        subplot(2,3,2);
        scatter(brevr,bievr);
        vline(1);
        xlim([-1.5,1.5]);
        xlabel('Real');
        ylabel('Imaginary')
        title('Eigenspectrum');
        axis equal
        box off
        grid()
        
        subplot(2,3,4:6);
        plot(Rates);
        grid();
        xlabel('Time (ms)');
        ylabel('Rate (Hz)');
        xlim([0,10000]);
        ylim([0,50]);
        box off
    
        subplot(2,3,1);
        imagesc(W);
        title('Connectivity matrix - Sparsity =',num2str(Sparsity));
        box off
    
        subplot(2,3,3);
        scatter(scores(:,1),scores(:,2));
        xlabel('PC1');
        ylabel('PC2');
        xlim([-300,300]);
        ylim([-300,300]);
        grid()
        box off
    
        set(gcf, 'WindowState', 'maximized');
    
        drawnow;  % Ensure the plot is fully rendered
        frame = getframe(gcf);  % Get the current frame
        im = frame2im(frame); % Convert the frame to an image
        imwrite(im, 'dynamics.png');

    else
        bevr = eig(W);
        [brevr,bixev] = sort(real(bevr),'descend');
        bievr = imag(bevr(bixev));
        [val,ix] = max(real(bevr));
        DominantMode = bevr(ix);


        figure;
        scatter(brevr,bievr);
        vline(1);
        xlim([-1.5,1.5]);
        xlabel('Real');
        ylabel('Imaginary')
        title('Eigenspectrum');
        axis equal
        box off

        set(gcf, 'WindowState', 'maximized');
        drawnow;  % Ensure the plot is fully rendered
        frame = getframe(gcf);  % Get the current frame
        im = frame2im(frame); % Convert the frame to an image
        imwrite(im, 'Spectrum.png');
        close all

    end
end
  
