function PCAPlot(Rates)

    [~,coeffs] = pca(Rates);
    plot(coeffs(:,1),coeffs(:,2));
    box off
    axis equal



end

