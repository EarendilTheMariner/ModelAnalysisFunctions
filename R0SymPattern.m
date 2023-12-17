function Stats = R0SymPattern(varargin)

    SLSweep = linspace(3000,50000,80);


     for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'SegmentLengthSweep'
                SLSweep = varargin{ii+1};
        end

    HeatMap = zeros(length(SLSweep),689);

    end
    for i = 1:length(SLSweep)
        close all
        clearvars -except SLSweep i HeatMap

        NetParams = NetworkParameters();
        NetParams.AddCelltype('V2a-1','bi','','','','','',1000,0.1*1.8);
        NetParams.AddCelltype('V1','bi','','','','','',3300,0.1*6);
        MouseNet = Network(SLSweep(i),1,NetParams);
        MouseNet.Simulate(8000);
        GeneratePlot(MouseNet);

        SSRates = MouseNet.Rates(end,:);
        [~,idx] = sort(MouseNet.Position(:,2));
        SSRates = SSRates(idx);
        
        HeatMap(i,:) = SSRates;

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
        bar(MouseNet.Position(:,2)/SLSweep(i),MouseNet.Rates(end,:),1,'histc');
    %    vline(SLSweep(i),'r','Segment boundary')
        xlabel('Norm. Rostro-caudal coordinate');
        ylabel('Activity (Firing rate)');
        xlim([0,1]);
        ylim([0,55]);
        legend('Segment Length =',num2str(round(SLSweep(i))));
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

        filename = 'R0SegmentLengthSweep.gif'

        % Write to GIF
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.4);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.4);
        end
    end
    save('Stats',"HeatMap");

    delete('temp.png');
end
