function EigenSpectrum(W,varargin)

  %  fig = figure;
    Eigenvalues = eig(W);
    [brevr,bixev] = sort(real(Eigenvalues),'descend');
    bievr = imag(Eigenvalues(bixev));
 %   fig = figure;
    scatter(brevr,bievr,'filled','MarkerFaceColor',Colors().BergBlack,'MarkerFaceAlpha',0.5)
    xlabel('Real','FontSize',20);
    ylabel('Imaginary','FontSize',20);
    vline(1);
  %  xlim([-0.3,1.3]);
  %  ylim([-1.3, 1.3]);
   
    axis equal
    box off
    set(gcf, 'WindowState', 'maximized');
    drawnow;  % Ensure the plot is fully rendered
%    save2pdf(fig,['./'],['Spectrum'],'-dpdf');
        

end

