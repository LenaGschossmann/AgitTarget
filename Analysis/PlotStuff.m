






%% Plot-Stuff
figure, plot(nS_pTrial, 'b')
figure, plot(meanAmp_pTrial, 'g')

figure, plot(1:length(trials_MSextract(itrial).left.Microsaccades.Amplitude(:)), trials_MSextract(itrial).left.Microsaccades.Amplitude(:))

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
