%calculate velocities
function [velosx velosy] = FEMvelo(trials_MSextract, lambda, sampleRate, degrees)

for itrial = 1:ntrials_tot
    
    lambda = 4; %after Englbert&Kliegl, 2002, but arbitrary?
    sampleRate = 0.250;
    degrees = 1;
    
    velosX = [];
    velosY = [];
    ind = 1;
    
    samples_x = trials_MSextract(itrial).left.samples.x(logical(trials_MSextract(itrial).left.samples.Good_Values(:)));
    samples_y = trials_MSextract(itrial).left.samples.y(logical(trials_MSextract(itrial).left.samples.Good_Values(:)));
    
    if degrees
        % distance from center in degrees (center is 0 degree)
        samples_x = -(1/experimentmat.px_per_deg).*(0.5*experimentmat.window_width - samples_x);
        samples_y = (1/experimentmat.px_per_deg).*(0.5*experimentmat.window_height - samples_y);
    end
    
    for ii = 3:(length(samples_x)-2) %because moving window -/+2
        velosX(ind) = (sum(samples_x(ii+1:ii+2))-sum(samples_x(ii-2:ii-1))) / (6*sampleRate); %Why 6, sampleRate = 250
        velosY(ind) = (sum(samples_y(ii+1:ii+2))-sum(samples_y(ii-2:ii-1))) / (6*sampleRate);
        ind = ind+1;
    end
    
    %median estimator (?) vgl Englbert&Kliegl 2002
    median_est_x = mean((velosX).^2) - (mean(velosX))^2; % ??
    std_x = sqrt(sum((velosX-median_est_x).^2)/(length(velosX)-1)); % ??
%     component_x = lambda*median_est_x;
    
    median_est_y = mean((velosY).^2) -(mean(velosY))^2; % ??
    std_y = sqrt(sum((velosY-median_est_y).^2)/(length(velosY)-1)); % ??
%     component_y = lambda*median_est_y;

    

    figure
    plot(velosX,velosY)
    t=-pi:0.01:pi;
    rx = std_x*lambda;
    ry = std_y*lambda;
    x=median_est_x+rx*cos(t);
    y=median_est_y+ry*sin(t);
    hold all, plot(x,y, 'black')    
end


end



%samples_no blinks


figure
g=gramm('x', velosX, 'y', velosY)
g.geom_line();
g.draw()


