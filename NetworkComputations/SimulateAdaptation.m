function Rates = SimulateAdaptation(W,t_steps,varargin)
    f_V_func = @tanh_f_V;
    tau_V = ones(size(W,1),t_steps)*50;% Slower Excitation
    %tau_V(sum(W,1) < 0,:) = 20; % Faster Inhibition
    noise_ampl = 0.4;
    seed = 5;
    I_e = ones(size(W,1),t_steps)*20;
    V_init = 0;
    threshold = 20;
    gain = ones(size(W,1),t_steps).*1;
    gainfeed = 0.2;
    f_max = 50;
    dfeed = 15;
    Source = ones(size(W,1));
    Target = zeros(size(W,1));
    Position = ones(size(W,1));
    Save = 0;
    dt = 10;
    A = zeros(size(W,1),size(W,2));
    D = zeros(size(W,1),size(W,1));

    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'tau_V'
                tau_V = varargin{ii+1};
            case 'noise_ampl'
                noise_ampl = varargin{ii+1};
            case 'seed'
                seed = varargin{ii+1};
            case 'I_e'
                
                I_e = varargin{ii+1};
            case 'V_init'
                V_init = varargin{ii+1};
            case 'threshold'
                threshold = varargin{ii+1};
            case 'gain'
                gain = varargin{ii+1};
            case 'fmax'
                f_max = varargin{ii+1};
            case 'Position'
                Position = varargin{ii+1};
            case 'ClassSource'
                Source = varargin{ii+1};
                LeftSource = Source & Position;
                RightSource = Source & ~Position;
            case 'ClassTarget'
                Target = varargin{ii+1};
                LeftTarget = Target & Position;
                RightTarget = Target & ~Position;
            case 'gainfeed'
                gainfeed = varargin{ii+1};    
        end
    end

    rng(seed);
    N = size(W,1);
    R = zeros(t_steps, N);  
    a = zeros(N, t_steps);
    V = zeros(N, t_steps);
    V(:, 1) = normrnd(V_init,3,size(V(:, 1)));
    R(1, :) = zeros(size(V_init));
    I_noise = normrnd(0.,noise_ampl,[1 (t_steps)*N]);
    I_rec = zeros(N, t_steps);
    I_tot = zeros(N, t_steps);

%    MR = (obj.Types == 'NF3')&(obj.Latera>0);
%    FlexR = (obj.MnID == 'Iliopsoas')&(obj.Latera>0);
%    ML = (obj.Types == 'NF3')&(obj.Latera<0);
%    FlexL = (obj.MnID == 'Iliopsoas')&(obj.Latera<0);

   
    SR = size(R);
    ColIx = [1:size(R,2)];
    tic
    for t = 2:t_steps
       for i = 1:N
            Rix = max(t-D(i,:)-1,1) + (ColIx-1)*SR(1);
            I_rec(i,t) = dot(W(i,:)+A(i,:), R(Rix));
            I_tot(i,t) = I_rec(i,t) + I_noise(t*i) + I_e(i, t);
            dV = (1. / tau_V(i, t)) * (-V(i,t-1) + I_tot(i,t));
            da = (1./ tau_a(i,t))*(I_rec(i,t)./nnz(W(i,:)));
            V(i,t) = V(i,t-1) + dV;
            a(i,t) = a(i,t-1) + da;
            if(Ad)
                A(i,:) = a(i,t)*W(i,:);
            end

            R(t, i) = f_V_func(V(i,t), threshold(i,t), gain(i, t), f_max(i));
            end
       end 

    
    toc
    Rates = R;
end

%% Helper Functions 
function val = tanh_f_V(V,threshold,gain,fmax) %20,1,100
    if (V-threshold)<=0.
        f = 0;
        %f = threshold*tanh(gain*(V-threshold)./threshold);
    elseif (V-threshold)>0
        f = fmax*tanh(gain*(V-threshold)./fmax);      
    end
    val = f;%+threshold;
end
function val = sigm_f_V(V,threshold,gain,fmax)
    val = fmax/(1+exp(-gain*(V-threshold)));
end
% function val = sigm_f_V(V,threshold,gain,fmax)
%     val = fmax/(1+exp(-3.9*gain*(V-threshold)/fmax));
% end


