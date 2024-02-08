function [WCtrl,DiffCtrl,EigenvaluesCtrl,RatesCtrl] = BumpControl(W,Coords,I,E,BumpIx,G)

    
    I = logical(I);
    E = logical(E);
    Dist = bsxfun(@minus,Coords,Coords');
    L = range(Coords,2);
    edges = -L:10:L;
    
    % Synapses belonging to inhibitory negative bump
    
    InhNegBump = (Dist(:,:) < 0) & (W(:,:) < 0);
    
    % Synapes belonging to inhibitory positive bump
    
    InhPosBump = (Dist(:,:) > 0) & (W(:,:) < 0); 
    
    % Synapse belonging to excitatory negative bump
    
    ExNegBump = (Dist(:,:) < edges(51)) & (W(:,:) > 0);
    
    % Synapses belonging to excitatory positive bump
    
    ExPosBump = (Dist(:,:) > edges(151)) & (W(:,:) > 0);
    
    % Synapses belonging to excitatory central bump
    
    CenterBump = (Dist(:,:) >= edges(50)) & (Dist(:,:) <= edges(150)) & (W(:,:) > 0);

    switch BumpIx
        case 1
            BIx = ExNegBump;
        case 2
            BIx = InhNegBump;
        case 3
            BIx = CenterBump;
        case 4
            BIx = InhPosBump;
        case 5
            BIx = ExPosBump;
        otherwise
            error('Invalid BumpIndex. Please choose a value between 1 and 5.');
    end

    % Outputs
    
    W(BIx) = G.*W(BIx);

    WCtrl = W;

%    WCtrl = BalanceConnectivity(WCtrl);

    WCtrl = BalanceNormalize(WCtrl);

    EigenvaluesCtrl = eig(WCtrl);
    [~,eix] = max(real(EigenvaluesCtrl));
    
    [brevr,bixev] = sort(real(EigenvaluesCtrl),'descend');
    bievr = imag(EigenvaluesCtrl(bixev));


    if real(EigenvaluesCtrl(eix)) >= 1

        RatesCtrl = SimulateNetwork(WCtrl,20000);
        RatesCtrl = RatesCtrl(10001:end,:);

    else
        RatesCtrl = [];
    end

    DiffCtrl = SpatialCoupling(WCtrl,E,I,Coords);


    % Visual overview
    
    fig = figure();
    
    subplot(2,2,1);
    bar(DiffCtrl);
    ylim([-30 70]);
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
    ylim([0,50]);
    box off
    
    set(gcf, 'WindowState', 'maximized');
    drawnow;  % Ensure the plot is fully rendered
    save2pdf(fig,['./'],['CtrlSummary'],'-dpdf');

end



