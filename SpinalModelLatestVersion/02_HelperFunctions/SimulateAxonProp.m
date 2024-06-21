function N = SimulateAxonProp(N,Type)
      %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Axon Growth Simulation Parameters
    it = 5000; % Number of iterations of simulations
    GrowthRate = 200;
    PosNoise = 0.1; %sqrt(DiffusionKernelWidthSR); % Noise on position (randon walk)

    % Video Related Parameters
    DisplayVid = 'on';
    isVideo = 0;
    mkdir(['./Videos']);
    v = VideoWriter(['./Videos' filesep 'AxonsProp_' PopOfInterest '-LR']);
    v.FrameRate = 30;
    v.Quality = 50;
    if(isVideo)
        open(v);
        DisplayVid = 'off';
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Select Points of pop with higher diff expression in space %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define region of interest
    ZeroMap = bwmorph(ZeroMap,'shrink',2);
    ZeroMap = bwmorph(ZeroMap,'spur');
    Cord = find(ZeroMap>0);
    
    % Populate postion Matrices for further computations
    posi = Coords(true,:);
    pos = posi + normrnd(0,1,size(posi));
    pos = fliplr(pos);
    pos =  repmat(pos,[2,1]);
    true = repmat(true,[2,1]);

    % Populate inital trajectory 
    Iin =  normrnd(0,PosNoise,size(pos));

    [LSU,LSD] = bounds(Landscape(Cord));
    [FX,FY] = gradient(Landscape,2);
    gamma = 1/LSD;

    lb = (pos(:,1) >= 1);
    rb = (pos(:,1) <= size(Landscape,1));
    bb = (pos(:,2) >= 1);
    ub = (pos(:,2) <= size(Landscape,2));
    bound = lb&rb&bb&ub;
    pos(bound==0,:) = [];
    true(bound==0) = [];
    s = [];
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Run Simulation and Plot Trajectories %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure('Color',Colors().BergGray09,'Position',[0,0,screensize(3)/2, screensize(4)]);
    set(fig, 'Visible',DisplayVid);
    ax1 = axes();
    imagesc(ax1,Landscape,"AlphaData",1);
    colormap(ax1,'gray');
    set(gca, 'YDir','reverse')
    set(gca, 'Color','none')
    axis off
    clim([LSU LSD]);
    
    ax2 = axes();
    scatter(ax2,CRD(:,1),CRD(:,2),150,Col,'filled','MarkerFaceAlpha',0.1)
    hold on
    scatter(ax2,pos(:,2),pos(:,1),50*(1+normalize(wn(true,ip),'range')),'white','filled','MarkerFaceAlpha',0.01)
    colormap(ax2,'gray');
    set(gca, 'YDir','reverse')
    set(gca, 'Color','none')
    axis off
    clim([LSU LSD]);


    linkaxes([ax1,ax2])
    PosIx = 1:size(pos,1);
    for tt = 1:it
        C = [];
        Z = [];
        for pp = 1:size(pos,1)
            Coord = floor(pos(pp,:));
            if(tt > 1)
                I = Traj(pp,:,tt-1,popind)-pos(pp,:);
            else 
                I = Iin(pp,:);
            end
            C = [C;Coord];
            Z = [Z;Landscape(Coord(1),Coord(2))];
            LocGrad = [FY(Coord(1),Coord(2)),FX(Coord(1),Coord(2))]; 
            update = sum([pos(pp,:);(wn(true(pp),ip)*GrowthRate*gamma*(Z(pp)*(LSD-Z(pp))).*LocGrad) + I/8000 + normrnd(0,PosNoise,size(pos,2))],1); 
            pos(pp,:) = update;
        end

        lb = ZeroMap(sub2ind(size(ZeroMap),round(pos(:,1),0),round(pos(:,2),0)));
        bound = lb;
        pos(bound==0,:) = [];
        PosIx(bound==0) = [];
        Traj(PosIx,:,tt,popind) = pos; 
        s1 = scatter(C(:,2),C(:,1),10,'filled','MarkerFaceColor','white','MarkerFaceAlpha',0.4);
        drawnow;
        delete(s1);
        if(isVideo)
            f = getframe(fig);
            writeVideo(v,f);
        end
    end
    if(isVideo)
        close(v);
    end
    popind = popind +1;
end