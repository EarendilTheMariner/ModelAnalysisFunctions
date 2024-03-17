function Stats = R0SpatialPattern(varargin)

    SLSweep = linspace(3000,50000,80);
    LSSweep = linspace(10,0.05,50);
    GainSweep = linspace(1,10,50);


     for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'SegmentLengthSweep'
                SLSweep = varargin{ii+1};
            case 'LengthScaleSweep'
                LSSweep = varargin{ii+1};
            case 'GainSweep'
                GainSweep = varargin{ii+1};
        end
     end
        
    Stats = zeros(length(SLSweep),689);

    for i = 1:length(SLSweep)
        close all
        clearvars -except SLSweep i Stats LSSweep GainSweep

        NetParams = NetworkParameters();
        NetParams.AddCelltype('V2a-1','cau','','','','','',1000,0.17);
        NetParams.AddCelltype('V1','bi','','','','','',3300,0.17);
        MouseNet = Network(3300,1,NetParams);
        MouseNet.Simulate(10000);
        GeneratePlot(MouseNet);

        SSRates = MouseNet.Rates(end,:);
        [~,idx] = sort(MouseNet.Position(:,2));
        SSRates = SSRates(idx);
        
        Stats(i,:) = SSRates;

        close all

        uniqueFilename = ['MouseNet_',num2str(i), '.mat']; 

        save(uniqueFilename, 'MouseNet');
    
        figure;
        subplot(2,2,3:4);
        plot(MouseNet.Rates);
        xlabel('Time (ms)');
        ylabel('Firing rate');
        xlim([0,8000]);
        ylim([0, 50]);
        grid();
        box off
    
        subplot(2,2,1:2);
        bar(MouseNet.Position(:,2)/3300,MouseNet.Rates(end,:),1,'histc');
    %    vline(SLSweep(i),'r','Segment boundary')
        xlabel('Norm. Rostro-caudal coordinate');
        ylabel('Activity (Firing rate)');
        xlim([0,1]);
        ylim([0,55]);
        legend('Segment length =',num2str(round(SLSweep(i))));
        title('Steady-state spatial activity')
        grid();
        box off

        set(gcf, 'WindowState', 'maximized');

        drawnow;  % Ensure the plot is fully rendered
        frame = getframe(gcf);  % Get the current frame
        im = frame2im(frame); % Convert the frame to an image
        imwrite(im, 'temp.png');
        im = imread('temp.png');
        [imind,cm] = rgb2ind(im,256);  % Convert the image to indexed color
        
        filename = 'R0InputPatternSweep.gif';
        
        % Write to GIF
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.4);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.4);
        end
    end
    save('Stats',"Stats");

    delete('temp.png');
end
