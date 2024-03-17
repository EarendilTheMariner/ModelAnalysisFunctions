function [WCtrl,DiffCtrl,EigenvaluesCtrl,RatesCtrl] = BumpControl(W,Coords,I,E,BumpIx,G,varargin)

    Balance = true;

    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'Balance'
                Balance = varargin{ii+1};
        end
    end

    
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

    Center1 = (Dist(:,:) >= edges(50)) & (Dist(:,:) <= edges(100)) & (W(:,:) > 0);
    Center2 = (Dist(:,:) <= edges(150)) & (Dist(:,:) >= edges(101)) & (W(:,:) > 0);

    switch BumpIx
        case 1
            BIx = ExNegBump;
        case 2
            BIx = InhNegBump;
        case 3
            BIx = Center1;
        case 4
            BIx = CenterBump;              
        case 5
            BIx = Center2;
        case 6
            BIx = InhPosBump;
        case 7
            BIx = ExPosBump;
        otherwise
            error('Invalid BumpIndex. Please choose a value between 1 and 7.');
    end

    % Outputs

    B2 = sum(W(InhNegBump));
    B4 = sum(W(InhPosBump));
    S = B2 + B4;


    if isequal(BIx, InhNegBump)

        W(BIx) = G.*W(BIx);
        
        G2 = (S-G*B2)/B4;

        W(InhPosBump) = G2*W(InhPosBump);

    elseif isequal(BIx,InhPosBump)

        W(BIx) = G.*W(BIx);

        G2 = (S-G.*B4)/B2;

        W(InhNegBump) = G2.*W(InhNegBump);

    else
        W(BIx) = G*W(BIx);
    end

    WCtrl = W;

    if Balance
        WCtrl = BalanceConnectivity(WCtrl);
    else
        disp('Unbalanced')
    end

    EigenvaluesCtrl = eig(WCtrl);
    [~,eix] = max(real(EigenvaluesCtrl));
    


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
    SpatialCoupling(WCtrl,E,I,Coords);
    
    ylim([-30 70])
    box off
    
    subplot(2,2,2);
    EigenSpectrum(WCtrl);

    
    subplot(2,2,3:4);
    plot(RatesCtrl);
    xlabel('Time','FontSize',20);
    ylabel('Rate','FontSize',20);
    xlim([0,10000]);
    ylim([0,70]);
    box off
    
    set(gcf, 'WindowState', 'maximized');
    drawnow;  % Ensure the plot is fully rendered
    save2pdf(fig,['./'],['CtrlSummary'],'-dpdf');

end



