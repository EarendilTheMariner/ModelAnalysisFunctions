function W = SymmetricBumps(W,E,I,Coords)

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
    
    Center1 = (Dist(:,:) >= edges(50)) & (Dist(:,:) <= edges(100)) & (W(:,:) > 0);
    Center2 = (Dist(:,:) <= edges(150)) & (Dist(:,:) >= edges(101)) & (W(:,:) > 0);

    %% Balancing Bump 1 and Bump 5
    balance = sum(W(ExNegBump))-sum(W(ExPosBump));

    bp = sum(W(ExNegBump))-(balance/2);
    bn = sum(W(ExPosBump))+(balance/2);
    
    facp = abs(bp/sum(W(ExNegBump)));
    
    facn = abs(bn/sum(W(ExPosBump)));

    W(ExNegBump) = W(ExNegBump)*facp;
    W(ExPosBump) = W(ExPosBump)*facn;

    %% Balancing Bump 2 and Bump 4

    balance = sum(W(InhNegBump))-sum(W(InhPosBump));


    bp = sum(W(InhNegBump))-(balance/2);
    bn = sum(W(InhPosBump))+(balance/2);
    
    facp = abs(bp/sum(W(InhNegBump)));
    
    facn = abs(bn/sum(W(InhPosBump)));

    W(InhNegBump) = W(InhNegBump)*facp;
    W(InhPosBump) = W(InhPosBump)*facn;

    %% Balance Bump 3

    balance = sum(W(Center1))-sum(W(Center2));

    bp = sum(W(Center1))-(balance/2);
    bn = sum(W(Center2))+(balance/2);

    facp = abs(bp/sum(W(Center1)));
    facn = abs(bn/sum(W(Center2)));

    W(Center1) = W(Center1)*facp;
    W(Center2) = W(Center2)*facn;


end

