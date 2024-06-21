function PlotTransitionRandomStructured(N,varargin)
    screensize = get(groot,'ScreenSize');
    InitPos = rand(size(N.Position));
    EndPos = N.Position;
    Pos = InitPos;
    S = N.Rates+10;
    S(S==0) = nan;
    Speed = 0.00005;
    noise_ample = 0.0005;    
    VA =[-180,0];
    t1 = 10;
    t2 = 3000;
    Save = 0;
    Mpos = mean(N.Position,1);

    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'SavePath'
                SavePath = varargin{ii+1};
                Save = 1;
            case 'SaveName'
                Name = varargin{ii+1};
        end
    end

    if (Save)
          if ~exist(SavePath, 'dir')
            mkdir(SavePath)
          end 
    v = VideoWriter([SavePath '/ ' Name],"MPEG-4");
    v.FrameRate = 30;
    v.Quality = 60;
    open(v);
    end
    
    DorPop = unique(N.Types(N.Layers == 'DRG'|N.Layers == '1Sp'| N.Layers == '2Sp0' | N.Layers == '2SpI' | N.Layers == '3Sp' | N.Layers == '4Sp'),'sorted');
    VentrPop = unique(N.Types(N.Layers == '5Spm'|N.Layers == '5SpL'| N.Layers == '6SpM' | N.Layers == '6SL' | N.Layers == '7Sp' | N.Layers == '8Sp' | N.Layers == 'D' | N.Layers == '10Sp' | N.Layers == 'Ps9'),'sorted');
    MN = unique(N.MnID(~isundefined(N.MnID)));

    cmapdorsal = fliplr(autumn(length(DorPop)));  %# Creates a 6-by-3 set of colors from the HSV colormap   
    cmapventralex = flipud(cool(length(VentrPop))); 
    cmapventralin = cool(length(VentrPop));
    cmapmn = sky(length(MN)); 

    C = zeros(size(N.ConnMat,1),3);

    for T = 1:length(DorPop)
        whr = N.Types == DorPop(T);
        C(whr,:) = repmat(cmapdorsal(T,:),[nnz(whr),1]);
    end
    
    for T = 1:length(VentrPop)
        whr = N.Types == VentrPop(T);
        if(mean(N.Transmit(whr)) < 0)
            C(whr,:) = repmat(cmapventralin(T,:),[nnz(whr),1]);
        else 
            C(whr,:) = repmat(cmapventralex(T,:),[nnz(whr),1]);
        end
    end  

    for T = 1:length(MN)
        whr = N.MnID == MN(T);
        C(whr,:) = repmat(cmapmn(T,:),[nnz(whr),1]);
    end

    fig = figure(Color=Colors().BergGray09,Position=screensize);
    ax1 = axes();
    axis off equal tight
    box off
    view(VA(1),VA(2))
    hold on
    %% Switch from a random network to a spatial network
    for j = 1:3:length(N.Rates)
        f1 = scatter3(ax1,Pos(:,1),Pos(:,2),Pos(:,3),0.65.*S(j,:)',C,'filled','MarkerFaceAlpha',0.75);
        drawnow
        if(Save)
           f = getframe(fig);
           writeVideo(v,f);
        end
        %% Update Position 
        if(j>t1)
            Pos = Pos + (EndPos-Pos)*Speed;% + normrnd(0.,noise_ample,[size(Pos)]);
        end
        if(j>t2)
            view([VA(1)-(90/(length(N.Rates)-t2))*(j-t2) VA(2)+(90/(length(N.Rates)-t2))*(j-t2)])
        end
        delete(f1);
    end
    view([-270 90]);
    %% Bring All neurons to same size
    S1 = 0.5.*S(end,:)';
    for ll = 1:0.1:5
        Size(:) =  S1(:) - ll*(0.5.*S(end,:)'- 2)./5;
        f1 = scatter(ax1,N.Position(:,1),N.Position(:,2),Size(end,:),C,'filled','MarkerFaceAlpha',0.75);
        if(Save)
           f = getframe(fig);
           writeVideo(v,f);
        end
        if(ll<5)
            delete(f1);
        end
    end

    T2 = text(Mpos(:,1),Mpos(:,2),'Rostro-caudal Synapse Distributions','FontSize',25,'Color',[1 1 1],'FontWeight','bold',VerticalAlignment='bottom');
    %% Display Properties of individual populations
    dis = 1;
    for T = {'V1','V2b','V2a-1','V2a-2','V3','DI6','MN','Rorb-I','RoraAdarb2','NF4','NF5'}
        whr = (N.Types == T{1}) & N.Segment=='L4';
        Mp = mean(N.Position(whr,:),1);
        SizeT = Size; 
        SizeT(whr) = SizeT(whr) + 50;
        for w = 1:10:600
           if(Save)
                f = getframe(fig);
                writeVideo(v,f);
            end
        end
        if(dis)
            delete(T2);
            dis = 0;
        end
        f1 = scatter(ax1,N.Position(:,1),N.Position(:,2),SizeT,C,'filled','MarkerFaceAlpha',0.75); 
        hold on 
        T1 = text(Mp(:,1),Mp(:,2),T{1},'FontSize',40,'Color',[1 1 1],'FontWeight','bold');
        for w = 1:10:600
           if(Save)
                f = getframe(fig);
                writeVideo(v,f);
            end
        end        
        GenerateDistribPlot(N,whr,mean(C(whr,:),1))
        delete(T1);
    end
    if(Save)
       f = getframe(fig);
       writeVideo(v,f);
       close(v)
    end
end

function GenerateDistribPlot(N,whr,C)
    ConnMat = N.ConnMat;
    Proj = logical(sum(ConnMat(:,whr),2)); 
    L = round(N.Position(:,2),-2);
    edges = [min(L):100:max(L)];
    yin = discretize(N.Position(Proj,2), edges);
    [~,~,Ci] = grpstats(Proj(Proj~=0), yin,["gname","mean","numel"]);

    if(mean(N.Transmit(whr)>0))
        mul = -1;
    else
        mul = 1;
    end
    yin(isnan(yin)) = [];
    barh(edges(unique(yin)+1), mul*10*Ci,'FaceColor',C,'FaceAlpha',0.5);
end
