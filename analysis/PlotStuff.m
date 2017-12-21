






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
template = logical(trials_MSextract(itrial).left.samples.Blink_Indices(:));
figure, plot(1:length(trials_MSextract(itrial).left.samples.x(template)), trials_MSextract(itrial).left.samples.x(template))


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
        g = gramm('x',log10(data.Amplitude),'color',data.condition);
        g.stat_density()
        g.draw()
        
        %% number of MS per trial
        stat = grpstats(data,{'trial','condition'});
        figure
        g = gramm('x',stat.condition,'y',stat.GroupCount,'color',stat.condition);
        g.stat_violin()
        g.draw()
        
        
%         
%         figure
%         g = gramm('x', data.DeltaX(ismember(data.subject, '7')), 'y', data.DeltaY(ismember(data.subject, '7')));
%         g.stat_smooth()
%         g.draw()


