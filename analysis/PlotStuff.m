
%plot microsaccade number against trials
figure
plot([1:80], tTrial.numMS)
condCol = {'r', 'g', 'b', 'y'}; 
tblock = [1:10;11:20;21:30;31:40;41:50;51:60;61:70;71:80];
for ii = 1:length(experimentmat.condition_list)
    hold all, plot(tblock(ii,:), tTrial.numMS(tTrial.trial(tblock(ii,:))),condCol{experimentmat.condition_list(ii)});
end


%plot x-velo against y-velo
figure
plot(trials_MSextract(itrial).left.samples.x, trials_MSextract(itrial).left.samples.y)
for ii = 1:length(trials_MSextract(itrial).left.Microsaccades.Start)
    hold all, plot(trials_MSextract(itrial).left.samples.x(trials_MSextract(itrial).left.Microsaccades.Start(ii):trials_MSextract(itrial).left.Microsaccades.End(ii)), trials_MSextract(itrial).left.samples.y(trials_MSextract(itrial).left.Microsaccades.Start(ii): trials_MSextract(itrial).left.Microsaccades.End(ii)), 'r')
end


%% Plot-Stuff
figure, plot(nS_pTrial, 'b')
figure, plot(meanAmp_pTrial, 'g')

figure, plot(trials_MSextract(itrial).left.Microsaccades.Start(:), trials_MSextract(itrial).left.Microsaccades.Amplitude(:))

figure, plot(number_blinks)

figure, plot(area_dispEll)

%% trials_work
%plot x and y values of single trials (set itrial)
figure, plot(1:length(trials_work(itrial).left.samples.x(:)), trials_work(itrial).left.samples.x(:))
% hold all, plot(1:length(trials(itrial).left.samples.x(:)), trials(itrial).left.samples.y(:), 'r')


for itrial = 1:10
    figure, plot(trials_work(itrial).left.samples.x, trials_work(itrial).left.samples.y)
    hold all, plot(trials_work(itrial).left.samples.x(trials_MSextract(itrial).left.Microsaccades.Start), trials_work(itrial).left.samples.y(trials_MSextract(itrial).left.Microsaccades.Start), 'r')
    hold all, plot(min(trials_work(itrial).left.samples.x) + range(trials_work(itrial).left.samples.x)/2, min(trials_work(itrial).left.samples.y) + range(trials_work(itrial).left.samples.y)/2, 'k*')
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

for itrial = 1:10
    figure, plot(trials_MSextract(itrial).left.samples.x, trials_MSextract(itrial).left.samples.y)
    hold all, plot(trials_MSextract(itrial).left.samples.x(trials_MSextract(itrial).left.Microsaccades.Start), trials_MSextract(itrial).left.samples.y(trials_MSextract(itrial).left.Microsaccades.Start), 'r')
    hold all, plot(screenRes.width/2, screenRes.height/2, 'k*')
end


for itrial = 1:5
%noBlinks
figure, plot(1:length(trials_work(itrial).left.samples.x(:)), trials_work(itrial).left.samples.x(:))
% hold all, plot(1:length(trials_work(itrial).left.samples.x(:)), trials_work(itrial).left.samples.y(:))
end


%% Plot with gramm

        %% Temporary plotting tools
        
        data = SAC_data;
        
        data(data.Amplitude>10^5,:) = [];
        %% Pure sanity checks, I don't expect any differences here
        % Draw Main Sequence
        figure
        g = gramm('x',log10(data.Amplitude),'y',log10(data.vPeak),'color',data.condition);
        g.geom_point();
        g.facet_grid([],data.condition)
        g.draw()
        %%
        figure
        g = gramm('x',log10(data.Amplitude),'y',log10(data.vPeak),'color',data.condition);
        g.stat_smooth();
        g.draw()
        
        %% Now let's have a look at more interesting things
        % amplitude densities
        figure
        g = gramm('x',data.Amplitude,'color',data.condition);
        g.stat_density()
        g.draw()
        
        %ellipse
        figure
        g=gramm('x', data.DeltaX, 'y', data.DeltaY, 'color', data.condition);
        g.stat_ellipse();
        g.draw()
        
        
        
        
        
        
        
        %% number of MS per trial

%         data.condition(strcmp(data.condition, 'CrossedBulleye')) = {'1'};
%         data.condition(strcmp(data.condition, 'GaussCirclesMasked_Dot')) = {'2'};
%         data.condition(strcmp(data.condition, 'StatCircles')) = {'3'};
%         data.condition(strcmp(data.condition, 'Bulleye')) = {'4'};
%         data.condition = cell2mat(data.condition);
%         data.condition = str2num(data.condition);
%         data.subject = str2num(data.subject);
        
%         stat = grpstats(data,'condition');
        figure
        g = gramm('x',data.condition,'color',data.condition);
        g.set_names('color', 'Condition');
        g.stat_bin('nbins', 4, 'width', 3);
        g.draw()
        
        
        %%plot eyemovement across all trials
        figure
        g = gramm('x', data.DeltaX, 'y', data.DeltaY, 'color', data.condition)
        g.geom_point();
        g.geom_vline('xintercept', 1920);
        g.geom_hline('yintercept', 1080);
        g.facet_grid([], data.condition);
        g.draw()
        
        %boxplot
        figure
        g= gramm('x', data.condition, 'y', data.Amplitude, 'color', data.condition);
        g.stat_summary();
        g.draw();


