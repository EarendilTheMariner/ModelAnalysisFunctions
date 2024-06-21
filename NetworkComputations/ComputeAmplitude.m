function MeanRMS = ComputeAmplitude(Rates)
    MR = mean(Rates);
    MR = repmat(MR,length(Rates(1,:)));
    Rates = Rates-MR;
    RMS = rms(Rates);
    MeanRMS = mean(RMS);
end

