function EigenSpectrum(W)

    
    Eigenvalues = eig(W);
    [brevr,bixev] = sort(real(Eigenvalues),'descend');
    bievr = imag(Eigenvalues(bixev));
    scatter(brevr,bievr);
    vline(1);
    % xlim([-1.5,1.5]);
    axis equal
    box off
    set(gcf, 'WindowState', 'maximized');


end

