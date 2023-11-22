function [dominantFrequency, eigFreq, lengthRatio] = ComputeFrequency(network_obj)

    %% function to compute the oscillatory frequency of the principal rate trace from network simulation

    [coeffs, scores] = pca(network_obj.Rates);
    pc1_scores = scores(:,1);

    length_ = length(network_obj.Rates(:,1));
    powers = fft(pc1_scores);

    powerNorm = abs(powers/length_);
    
    pos_freq = powerNorm(1:(length_/2));

    pos_freq(2:end-1) = 2*pos_freq(2:end-1);

    dt = 0.001;

    nyquist = 1/(2*dt); % 500 cycles per second

    freqs = linspace(0, nyquist, length_/2);

    scoresOnesided = pos_freq;

    [~, idx] = max(scoresOnesided);

    dominantFrequency = freqs(idx);

    eigs = eig(network_obj.ConnMat);
    eigIx = real(eigs) > 1;
    eigFreq = imag(eigs(eigIx))/2*pi;

    lengthRatio = network_obj.Parameters.LengthScales(1)/network_obj.Parameters.LengthScales(2); % exc. to inh. lengthscale ratio

    
end


     


