
%plot microsaccade number against trials
figure
plot([1:80], tTrial.numMS)
condCol = {'r', 'g', 'b', 'y'}; 
tblock = [1:10;11:20;21:30;31:40;41:50;51:60;61:70;71:80];
for ii = 1:length(experimentmat.condition_list)
    hold all, plot(tblock(ii,:), tTrial.numMS(tTrial.trial(tblock(ii,:))),condCol{experimentmat.condition_list(ii)});
end


%plot x against y
figure
plot(trials_MSextract(itrial).left.samples.x, trials_MSextract(itrial).left.samples.y)
for ii = 1:length(trials_MSextract(itrial).left.Microsaccades.Start)
    hold all, plot(trials_MSextract(itrial).left.samples.x(trials_MSextract(itrial).left.Microsaccades.Start(ii):trials_MSextract(itrial).left.Microsaccades.End(ii)), trials_MSextract(itrial).left.samples.y(trials_MSextract(itrial).left.Microsaccades.Start(ii): trials_MSextract(itrial).left.Microsaccades.End(ii)), 'r')
%     t=-pi:0.01:pi;
%     rx = std_x*lambda;
%     ry = std_y*lambda;
%     x=median_est_x+rx*cos(t);
%     y=median_est_y+ry*sin(t);
%     hold all, plot(x,y, 'black')
end


%% Plot-Stuff
figure, plot(nS_pTrial, 'b')
figure, plot(meanAmp_pTrial, 'g')

figure, plot(trials_MSextract(itrial).left.Microsaccades.Start(:), trials_MSextract(itrial).left.Microsaccades.Amplitude(:))

figure, plot(number_blinks)

figure, plot(area_dispEll)

%% trials_raw
%plot x and y values of single trials (set itrial)
figure, plot(1:length(trials_raw(itrial).left.samples.x(:)), trials_raw(itrial).left.samples.x(:))
% hold all, plot(1:length(trials(itrial).left.samples.x(:)), trials(itrial).left.samples.y(:), 'r')

%%
for itrial = [10 11]
    figure, plot(trials_MSextract(itrial).left.samples.x, trials_MSextract(itrial).left.samples.y)
    hold all, 
    plot(trials_MSextract(itrial).left.samples.x(trials_MSextract(itrial).left.Microsaccades.Start), trials_MSextract(itrial).left.samples.y(trials_MSextract(itrial).left.Microsaccades.Start), 'r')
    plot(min(trials_MSextract(itrial).left.samples.x) + range(trials_MSextract(itrial).left.samples.x)/2, min(trials_MSextract(itrial).left.samples.y) + range(trials_MSextract(itrial).left.samples.y)/2, 'k*')
    plot(3840/2,2160/2,'k*')
end




%% MSextract
%plot x and y values of single trials (set itrial)
figure, plot(1:length(trials_MSextract(itrial).left.samples.x(:)), trials_MSextract(itrial).left.samples.x(:))
% hold all, plot(1:length(trials(itrial).left.samples.x(:)), trials(itrial).left.samples.y(:), 'r')

%samples_no blinks
template = logical(trials_MSextract(itrial).left.samples.Good_Values(:));

figure, plot(1:length(trials_MSextract(itrial).left.samples.x(template)), trials_MSextract(itrial).left.samples.x(template))

%plot eyemovement
figure
g = gramm('x', trials_MSextract(itrial).left.samples.x(template), 'y', trials_MSextract(itrial).left.samples.y(template))
g.geom_point();
g.geom_vline('xintercept', 1920);
g.geom_hline('yintercept', 1080);
g.draw()
%%
for itrial = 19
    samp = trials_MSextract(itrial).left.samples;
    micro = trials_MSextract(itrial).left.Microsaccades;
    
    % modified jet-colormap



    

    figure, 
    p = plot(samp.x, samp.y,'o-b');
    drawnow
    n = ceil(length(samp.x)/50);
    
    cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
    cd = repmat(cd,1,50)
    cd = cd(:,1:length(samp.x));
    
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    hold all, 
    line([samp.x(micro.Start);samp.x(micro.End)], [samp.y(micro.Start);samp.y(micro.End)],'Color','red','LineWidth',2)
    hold all, plot(screenRes.width/2, screenRes.height/2, 'k*')
end
%%

for itrial = 1:5
%noBlinks
figure, plot(1:length(trials_raw(itrial).left.samples.x(:)), trials_raw(itrial).left.samples.x(:))
% hold all, plot(1:length(trials_raw(itrial).left.samples.x(:)), trials_raw(itrial).left.samples.y(:))
end


%% Plot with gramm

        %% Temporary plotting tools
        
        data = SAC_data;
        data.subject = str2num(data.subject);
        if size(data, 1) == size(MS_data, 1)
            msORsac = 'Microsaccades after E&K';
        else
            msORsac = 'Saccades after Eyelink';
            data.DeltaX = abs(data.Sx - data.Ex);
            data.DeltaY = abs(data.Sy - data.Ey);
            if sum(data.DeltaX > 10e+05)>0
                data.DeltaX(data.DeltaX > 10e+05) = 0;
                data.DeltaY(data.DeltaX > 10e+05) = 0;
            end
            if sum(data.DeltaY > 10e+05)>0
                data.DeltaX(data.DeltaY > 10e+05) = 0;
                data.DeltaY(data.DeltaY > 10e+05) = 0;    
            end
            data.Amplitude = sqrt(data.DeltaX.^2 + data.DeltaY.^2);
            data.vPeak = data.Speed;
        endSAC

            
        %% Pure sanity checks, I don't expect any differences here
        % Draw Main Sequence
            figure
            g = gramm('x',(data.Amplitude),'y',log10(data.vPeak),'color',data.condition);
            g.set_title(sprintf('Microsaccades of Subjects %s', int2str(subjects)));
            g.geom_point('alpha',0.1);
            g.facet_grid(data.condition, data.subject, 'row_labels', false);
            g.axe_property('YLim', [1 2], 'XLim', [0 500])
            g.draw()
        
%         figure
%         g = gramm('x',log10(data.Amplitude),'y',log10(data.vPeak),'color',data.condition);
%         g.stat_smooth();
%         g.draw()
        

        %% Now let's have a look at more interesting things
        
        if size(data, 1) == size(MS_data, 1)
        % amplitude densities
            figure
            g = gramm('x',data.Amplitude(data.Amplitude<400),'color',data.condition(data.Amplitude<400));
            g.stat_density()
            g.facet_grid(data.subject(data.Amplitude<400),[])
            g.draw()
            %%
        %ellipse
        figure
        g=gramm('x', data.DeltaX,'y', data.DeltaY, 'color', data.condition);
        g.stat_ellipse();
        g.set_title(sprintf('Deviation from 0 [px] during Microsaccades of Subjects %s', int2str(subjects)));
        g.facet_grid([], data.subject)
        g.axe_property('Xlim',[-150 150],'Ylim',[-100 100],'DataAspectRatio',[1 1 1]);
        g.draw()
        end
        
        
        %% number of saccades

        figure
        g = gramm('x',data.condition,'color',data.condition);
        g.set_names('color', 'Condition');
        g.stat_bin('nbins', 4, 'width', 3);
        g.set_title(sprintf('# %s of Subjects %s', msORsac, int2str(subjects)));
        g.facet_grid(data.subject,[]);
        g.draw()
        
%         figure
%         g = gramm('x', data.Amplitude, 'color', data.condition);
%         g.set_names('color', 'Condition');
%         g.stat_bin('nbins', 50);
%         g.facet_grid(data.subject,[])
%         g.draw();
%         

        %% plot eyemovement across all trials
            figure
            g = gramm('x', data.DeltaX, 'y', data.DeltaY, 'color', data.condition)
            g.geom_point();
    %         g.geom_vline('xintercept', 1920);
    %         g.geom_hline('yintercept', 1080);
            g.facet_grid(data.condition, data.subject, 'row_labels', false);
            g.axe_property('Xlim',[-200 200],'Ylim',[-200 200],'DataAspectRatio',[1 1 1]);
            g.set_title(sprintf('X and Y dispersion of %s [px] of Subjects %s', msORsac, int2str(subjects)));
            g.draw()
       