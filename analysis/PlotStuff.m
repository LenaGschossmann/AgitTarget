
%plot microsaccade number against trials
figure
plot([1:80], tTrial.numMS)
condCol = {'r', 'g', 'b', 'y'}; 
tblock = [1:10;11:20;21:30;31:40;41:50;51:60;61:70;71:80];
for ii = 1:length(expmat.condition_list)
    hold all, plot(tblock(ii,:), tTrial.numMS(tTrial.trial(tblock(ii,:))),condCol{expmat.condition_list(ii)});
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




%% MSextract
%plot x and y values of single trials (set itrial)
figure, plot(1:length(trials_MSextract(itrial).left.samples.x(:)), trials_MSextract(itrial).left.samples.x(:))
hold all, plot(trials_MSextract(itrial).left.Microsaccades.Start(:), trials_MSextract(itrial).left.samples.x(trials_MSextract(itrial).left.Microsaccades.Start(:)), '*r');

%samples_no blinks
template = logical(trials_MSextract(itrial).left.samples.Good_Values(:));

figure, plot(1:length(trials_MSextract(itrial).left.samples.x(template)), trials_MSextract(itrial).left.samples.x(template))

for itrial = [10:18]
    template = logical(trials_MSextract(itrial).left.samples.Good_Values(:));
    figure, plot(trials_MSextract(itrial).left.samples.x(template), trials_MSextract(itrial).left.samples.y(template))
    hold all, 
    plot(trials_MSextract(itrial).left.samples.x(trials_MSextract(itrial).left.Microsaccades.Start), trials_MSextract(itrial).left.samples.y(trials_MSextract(itrial).left.Microsaccades.Start), 'r')
%     plot(min(trials_MSextract(itrial).left.samples.x) + range(trials_MSextract(itrial).left.samples.x)/2, min(trials_MSextract(itrial).left.samples.y) + range(trials_MSextract(itrial).left.samples.y)/2, 'k*')
    plot(3840/2,2160/2,'k*')
end


%plot eyemovement
figure
g = gramm('x', trials_MSextract(itrial).left.samples.x(template), 'y', trials_MSextract(itrial).left.samples.y(template))
g.geom_point();
g.geom_vline('xintercept', 1920);
g.geom_hline('yintercept', 1080);
g.draw()

%%Bene
for itrial = 19
    samp = trials_MSextract(itrial).left.samples;
    micro = trials_MSextract(itrial).left.Microsaccades;
    
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


for itrial = 1:5
%noBlinks
figure, plot(1:length(trials_raw(itrial).left.samples.x(:)), trials_raw(itrial).left.samples.x(:))
% hold all, plot(1:length(trials_raw(itrial).left.samples.x(:)), trials_raw(itrial).left.samples.y(:))
end




%% Plot with gramm
        
        data = MS_data;
        data.subject = str2num(data.subject);
        OV_conditions.subject = str2num(OV_conditions.subject);
        OV_trials.subject = str2num(OV_trials.subject);
        
%         label_sbj = [];
%         for iS = 1:length(subjects)
%             label_sbj = [label_sbj; sprintf('Subject %i', subjects(iS))];
%         end   
        
        if size(data, 1) == size(MS_data, 1)
            msORsac = 'Microsaccades detected after Engbert & Kliegl (2003)';
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
        end

            
        %% Pure sanity checks, I don't expect any differences here
        % Draw Main Sequence
            figure
            g = gramm('x',log10(data.Amplitude),'y',log10(data.vPeak),'color',data.condition);
            g.set_title(msORsac);
            g.geom_point('alpha',0.1);
            g.facet_grid(data.condition, data.subject, 'row_labels', false);
            g.set_names('column', 'Subject', 'x', 'Amplitude (log)', 'y', 'Peak velocity (log)');
%             g.axe_property('YLim', [1.5 2], 'XLim', [0 300]) %SAC
            g.axe_property('YLim', [2 4], 'XLim', [0 3]) %MS
            g.draw()
        
%         figure
%         g = gramm('x',log10(data.Amplitude),'y',log10(data.vPeak),'color',data.condition);
%         g.stat_smooth();
%         g.draw()
            
            % valid trials
            figure
            g = gramm('x', OV_conditions.valid_trials, 'color', OV_conditions.condition);
            g.stat_bin();
            g.facet_grid(OV_conditions.subject, []);
            g.set_names('x', 'Number of valid trials', 'row', 'Subject');
            g.draw()

        

        %% Now let's have a look at more interesting things    
        if size(data, 1) == size(MS_data, 1)
        % amplitude densities
            figure
            g = gramm('x',data.Amplitude,'color',data.condition);
            g.stat_density()
            g.axe_property('XLim', [0 75]);
            g.facet_grid(data.subject,[]);
            g.set_names('column', 'Subject', 'x', 'Amplitude [px]', 'row', 'Subject');
            g.set_title('Distribution of Microsaccade Amplitudes');
            g.draw()
            %ellipse
            figure
            g=gramm('x', (data.EndX-expmat.screenRes.width/2),'y', (data.EndY-expmat.screenRes.height/2), 'color', data.condition);
%             g.geom_point('alpha', 0.05);
            g.stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
            g.set_title('Deviation of Microsaccades from Fixation point');
            g.facet_grid([], data.subject);
            g.set_names('column', 'Subject', 'x', '[px]', 'y', '[px]');
            g.axe_property('Xlim',[-500 400],'Ylim',[-500 500],'DataAspectRatio',[1 1 1]);
            g.draw()        
        end
             
        %% number of saccades
%         figure
%         g = gramm('x',data.condition,'color',data.condition);
%         g.set_names('color', 'Condition');
%         g.stat_bin('nbins', 4, 'width', 3);
%         g.set_title(sprintf('# %s of Subjects %s', msORsac, int2str(subjects)));
%         g.facet_grid(data.subject,[]);
%         g.draw()
        
%         figure
%         g = gramm('x', data.Amplitude, 'color', data.condition);
%         g.set_names('color', 'Condition');
%         g.stat_bin('nbins', 50);
%         g.facet_grid(data.subject,[])
%         g.draw();
%         

        %% plot microsaccades across all trials
        figure
        g= gramm('x', [data.StartX;data.EndX], 'y',  [data.StartY;data.EndY], 'color', [data.condition; data.condition])
        g.geom_line('alpha', 0.08);
        g.geom_hline('yintercept', expmat.screenRes.height/2);
        g.geom_vline('xintercept', expmat.screenRes.width/2);
        g.facet_grid([data.condition; data.condition], [], 'row_labels', false);
        g.axe_property('DataAspectRatio',[1 1 1]);
        g.set_names('x', '[px]', 'y', '[px]');
        g.set_title(sprintf('%s across Subjects (crossing = Fixation point)', msORsac));
        g.draw()
        
        figure
        g= gramm('x', [data.StartX;data.EndX], 'y',  [data.StartY;data.EndY], 'color', [data.condition data.condition])
        g.geom_line('alpha', 0.08);
        g.geom_hline('yintercept', expmat.screenRes.height/2);
        g.geom_vline('xintercept', expmat.screenRes.width/2);
        g.facet_grid([data.condition; data.condition], [data.subject; data.subject], 'row_labels', false);
        g.set_names('column', 'Subject', 'x', '[px]', 'y', '[px]');
        g.axe_property('DataAspectRatio',[1 1 1]);
        g.set_title(sprintf('%s per Subject (crossing = Fixation point)', msORsac));
        g.draw()
      
            
        %% Number of MS
        
        OV_trials.BlockCond = repmat({'x'}, size(OV_trials, 1), 1); 
        
        blockID1 = strcmp(OV_trials.condition, 'GaussCirclesMasked_Dot')==1 & OV_trials.trial <= 40;
        blockID2 = strcmp(OV_trials.condition, 'CrossedBulleye')==1 & OV_trials.trial <= 40;
        blockID3 = strcmp(OV_trials.condition, 'Bulleye')==1 & OV_trials.trial <= 40;
        blockID4 = strcmp(OV_trials.condition, 'StatCircles')==1 & OV_trials.trial <= 40;
        blockID5 = strcmp(OV_trials.condition, 'GaussCirclesMasked_Dot')==1 & OV_trials.trial > 40;
        blockID6 = strcmp(OV_trials.condition, 'CrossedBulleye')==1 & OV_trials.trial > 40;
        blockID7 = strcmp(OV_trials.condition, 'Bulleye')==1 & OV_trials.trial > 40;
        blockID8 = strcmp(OV_trials.condition, 'StatCircles')==1 & OV_trials.trial > 40;
        
        OV_trials.BlockCond(blockID1) = repmat({'GaussCirclesMasked_Dot_1'},70,1);
        OV_trials.BlockCond(blockID2) = repmat({'CrossedBulleye_1'},70,1);
        OV_trials.BlockCond(blockID3) = repmat({'Bulleye_1'},70,1);
        OV_trials.BlockCond(blockID4) = repmat({'StatCircles_1'},70,1);
        OV_trials.BlockCond(blockID5) = repmat({'GaussCirclesMasked_Dot_2'},70,1);
        OV_trials.BlockCond(blockID6) = repmat({'CrossedBulleye_2'},70,1);
        OV_trials.BlockCond(blockID7) = repmat({'Bulleye_2'},70,1);
        OV_trials.BlockCond(blockID8) = repmat({'StatCircles_2'},70,1);
        
        figure
        g = gramm('x', OV_trials.trial , 'y', OV_trials.norm_numberMS, 'color', OV_trials.BlockCond, 'group', OV_trials.condition)
        g.stat_smooth();
        g.facet_grid(OV_trials.subject,[], 'scale','free_X');
        g.set_color_options('n_color', 8, 'n_lightness', 2, 'map', [0.9 0.4 0.4; 0 0 0; 0.9 0.4 0.4; 0 0 0; 0.5 0.7 0.1; 0 0 0; 0.5 0.7 0.1; 0 0 0; 0.2 0.8 0.9; 0 0 0; 0.2 0.8 0.9; 0 0 0; 0.4 0.3 0.7; 0 0 0; 0.4 0.3 0.7; 0 0 0]);
%         g.axe_property('DataAspectRatio', [1 1 1]);
        g.set_names('row', 'Subject');
        g.set_title(sprintf('Number of %s', msORsac));
        g.draw()
        

            
       