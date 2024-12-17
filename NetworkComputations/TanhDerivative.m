function val_deriv = TanhDerivative(V, threshold, gain, fmax)
    % Initialize f_deriv
    f_deriv = zeros(size(gain));
    
    % Calculate the derivative based on the condition
    neg = (V - threshold) <= 0;
    pos =  (V - threshold) > 0;
    % Derivative when V - threshold <= 0
    f_deriv(neg) = gain(neg) .* (1 - tanh(gain(neg) .* (V(neg) - threshold(neg))./ threshold(neg)).^2); %./ threshold(neg))
    % Derivative when V - threshold > 0
    f_deriv(pos) = gain(pos) .* (1 - tanh(gain(pos) .* (V(pos) - threshold(pos))./ fmax(pos)).^2); %./ fmax(pos))
    % The derivative of val with respect to V
    val_deriv = f_deriv;
end