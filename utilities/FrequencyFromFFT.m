function Frequency = FrequencyFromFFT(Signal)

    Fs = 1;                 % Sampling frequency                    
    T = 1/Fs;               % Sampling period       
    L = length(Signal);     % Length of signal
    t = (0:L-1)*T;          % Time vector

    Y = fft(Signal);
    X = Fs/L*(0:L-1);
    Y(length(X)/2:end) = 0;
    Y(1) = 0;


    [~,ix] = max(abs(Y));

    Frequency = X(ix);
    %Frequency = Frequency*1000;
    

  %  figure;
  %  plot(X,abs(Y),"LineWidth",3)
  %  title("Complex Magnitude of fft Spectrum")
  %  xlabel("f (Hz)")
  %  ylabel("|fft(X)|")
end


