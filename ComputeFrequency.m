function [DominantMode, Powers] = ComputeFrequency(Sig)

    %% 

    Length_ = length(Sig);
    Powers = fft(Sig);
    PowerNorm = abs(Powers/Length_);
    PosFreq = PowerNorm(1:(Length_/2));
    Nyquist = 0.5*Length_; % 1 cycle / 2 data point
    Freqs = linspace(0, Nyquist, Length_/2);
    [~, idx] = max(real(PosFreq));
    DominantMode = Freqs(idx);

    plot(Freqs, PosFreq)
    title('Magnitude Spectrum')
    xlabel('Normalized Frequency, k')
    ylabel('|F(k)|')
    
end


     


