
function [W,P] = spinal_network(size,f_exc,varargin)
    
    mean = 1;
    std = 0.1;

    % Parsing variable arguments 
    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'Mean'
                mean = varargin{ii+1};
            case 'Std'
                std = varargin{ii+1};
        end
    end

    n = size;
    W = zeros(n);
    exc_bias = 0.3;

    neuron_types = [ones(1, n*f_exc), -1*ones(1,(n*(1-f_exc)))]; % 1xn array signifying neuron types (1 or -1);
    caud_coords = 10 * rand(n, 1); %caudal coordinate
    lat_coords = 3 * rand(n,1); %lateral coordinate
    dors_coords = 3 * rand(n,1); %dorsal coordinate

    % pair wise distances for each coordinate (Euclidean)

    dist_caud = pdist2(caud_coords, caud_coords);
    dist_lat = pdist2(lat_coords, lat_coords);
    dist_dors = pdist2(dors_coords, dors_coords);

    %connection probabilities: exponential decay from coordinate distances

    connection_probs = exp(-((dist_caud+dist_lat+dist_dors)/3));

    % Weights for inhibitory projections
    for i = n*(1-f_exc)+1:n
        for j = 1:n
            if i ~= j && rand <= connection_probs(i,j)
                W(i,j) = -1*abs(normrnd(mean,std));

            end
        end
    end
    
    % Excitatory biases in the caudal directon

    for i = 1:n*f_exc
        for j = 1:n
            if caud_coords(j) > caud_coords(i)
                connection_probs(i,j) = connection_probs(i,j) + exc_bias;
            end
            if i ~= j && rand <= connection_probs(i,j)
                W(i,j) = +1*abs(normrnd(mean,std));
            end
        end
    end
    
   W = BalanceConnectivity(W'); %balance the matrix
   P = caud_coords;
end

  