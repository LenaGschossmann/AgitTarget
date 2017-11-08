
addpath('/net/home/student/l/lgschossmann/AgitTarget/edfread/build/linux64');

%Indicate current file
curSub = 'AgitT_sb_3';

%% .mat file

load(sprintf('/net/home/student/l/lgschossmann/AgitTarget/Data/%s.mat', curSub));

%% EDF data

[trials, info] = edfread(sprintf('/net/home/student/l/lgschossmann/AgitTarget/Data/%s.EDF', curSub), 'write');


%% Create working trials structure

trials_work = trials;
for itrial = 1:ntrials_tot
    if (isstruct(trials_work(itrial).left)), else, trials_work(itrial).left = trials_work(itrial).right; end
end


%% Cut off first and last 2 seconds (equals 250 samples) of each trial
cutoff_sec = 3;
cutoff = cutoff_sec*250;

for itrial = 1:ntrials_tot
    borders = [(cutoff+1) length(trials_work(itrial).left.samples.x)-cutoff];
    trials_work(itrial).left.samples.time = trials_work(itrial).left.samples.time(borders(1):borders(2));
    trials_work(itrial).left.samples.x = trials_work(itrial).left.samples.x(borders(1):borders(2));
    trials_work(itrial).left.samples.y = trials_work(itrial).left.samples.y(borders(1):borders(2));
    trials_work(itrial).left.samples.pupil = trials_work(itrial).left.samples.pupil(borders(1):borders(2));
end


%% cut out blinks
%blink detection: difference between two adjacent points extremely big

% number_blinks = zeros(1, ntrials_tot);
% 
% for itrial = 1:ntrials_tot
%     samples_x = trials_work(itrial).left.samples.x(:); samples_y =
%     trials_work(itrial).left.samples.y(:); pupil =
%     trials_work(itrial).left.samples.pupil(:);
%     
%     if(sum(samples_x(:) >= 1e07)>0) %check if there are extremely high
%     values in vector which indicate blinks
%         diffs_x = NaN(1, length(samples_x)-1);
% 
%         diffs_x(:) = abs(diff(samples_x));
% 
%         [diffs_sorted, sortIdx_x] = sort(diffs_x, 'descend');
% 
%         %number of blinks [M, blinks] = max(abs(diff(diffs_sorted)));
% 
%         blinksIdx = sort(sortIdx_x(1:blinks)); %indices of blinks in
%         samples_x in right chronological order
%         
%         %cut out area around blinks for blIdx = 1:blinks
%             if (mod(blIdx,2)>0)
%                 goOn = 1; while(goOn == 1)
%                     if ((samples_x(blinksIdx(blIdx)) -
%                     samples_x(blinksIdx(blIdx)-1)) >= 2)
%                         blinksIdx(blIdx) = blinksIdx(blIdx)-1;
%                     else
%                         goOn = 0;
%                     end
%                 end
%             else
%                 goOn = 1; while(goOn == 1)
%                     if ((samples_x(blinksIdx(blIdx))) -
%                     samples_x(blinksIdx(blIdx)+1) >= 2)
%                         blinksIdx(blIdx) = blinksIdx(blIdx)+1;
%                     else
%                         goOn = 0;
%                     end
%                 end
%             end
%         end
%         
%         %Replace blink values with interpolation if (blinks >1)
%             blIdx = 2; %Index of end of each blink while (blIdx <=
%             blinks)
%                 tmpVec = blinksIdx(blIdx-1)+1:blinksIdx(blIdx);
%                 samples_x(tmpVec) =
%                 linspace(samples_x(blinksIdx(blIdx-1)),
%                 samples_x(blinksIdx(blIdx)+1), length(tmpVec));
%                 samples_y(tmpVec) =
%                 linspace(samples_y(blinksIdx(blIdx-1)),
%                 samples_y(blinksIdx(blIdx)+1), length(tmpVec));
% 
% %                 % Replace blink-whole with value of last sample before
% %                 samples_x(tmpVec) = samples_x(blinksIdx(blIdx-1)); %
% samples_y(tmpVec)40 = samples_y(blinksIdx(blIdx-1));
%         %         blinksIdx(blIdx-1) = blinksIdx(blIdx-1)+1;  %blink
%         starts with first point with much higher value (and with
%         diff-function the index leads to the samplepoint before) %
%         blinksIdx_y(blIdx-1) = blinksIdx_y(blIdx-1)+1; %
%         blinksIdx_p(blIdx-1) = blinksIdx_p(blIdx-1)+1;
%  
%                 blIdx = blIdx+2;
%             end
%             
%             %check if there is any single blink at the end if
%             (mod(blinks, 2)>0)
%                 blIdx = blIdx-1; tmpVec =
%                 blinksIdx(blIdx)+1:length(samples_x); samples_x(tmpVec) =
%                 samples_x(blinksIdx(blIdx-1)); samples_y(tmpVec) =
%                 samples_y(blinksIdx(blIdx-1));
%             end
%             
%         else
%             blIdx = 1; tmpVec = blinksIdx(blIdx)+1:length(samples_x);
%             samples_x(tmpVec) = samples_x(blinksIdx(blIdx-1));
%             samples_y(tmpVec) = samples_y(blinksIdx(blIdx-1));
%             
%         end
% 
% 
%         trials_work(itrial).left.samples.x(:) = samples_x(:);
%         trials_work(itrial).left.samples.y(:) = samples_y(:);
%         trials_work(itrial).left.samples.pupil(:) = pupil(:);
%         number_blinks(itrial) = round(blinks/2);
%     end
%     
% end
% 
% clear('samples_x', 'samples_y', 'pupil', 'sortIdx_x', 'sortIdx_y',
% 'sortIdx_p', 'diffs_sorted', 'blIdx', 'blinksIdx', 'blinksIdx_y',
% 'blinksIdx_p', 'M');




%% Extract Microsaccades

[trials_MSextract]= edfExtractMicrosaccades(trials_work);
    % MS parameters
    %     --> Amplitude: in px
    %     --> Phi: Amplitude in vis deg
    %     --> normgx & normgy: distance in px
    %     --> vPeak: Peak velocity


%% Disregard all MS during blinks if blinks were not cut out before

number_blinks = zeros(1, ntrials_tot);

for itrial = 1:ntrials_tot
    Blink_Indices = ones(1, length(trials_work(itrial).left.samples.x)); %Contains 0 at every element that has to be disregarded as it lies in a blink
    
    samples_x = trials_work(itrial).left.samples.x(:);
    samples_y = trials_work(itrial).left.samples.y(:);
    pupil = trials_work(itrial).left.samples.pupil(:);
    
    if(sum(samples_x(:) >= 1e07)>0) %check if there are extremely high values in vector which indicate blinks
        diffs_x = NaN(1, length(samples_x)-1);

        diffs_x(:) = abs(diff(samples_x));

        [diffs_sorted, sortIdx_x] = sort(diffs_x, 'descend');

        %number of blinks
        [M, blinks] = max(abs(diff(diffs_sorted))); %blinks/2 = number of blinks
        blinksIdx = sort(sortIdx_x(1:blinks)); %indices of blinks in samples_x in right chronological order
        
        %include area around blinks
        for blIdx = 1:blinks
            if (mod(blIdx,2)>0)
                goOn = 1;
                while(goOn == 1)
                    if ((samples_x(blinksIdx(blIdx)) - samples_x(blinksIdx(blIdx)-1)) >= 2) % 2 is randomly chosen
                        blinksIdx(blIdx) = blinksIdx(blIdx)-1;
                    else
                        while (samples_x(blinksIdx(blIdx)-1) < 0)
                            blinksIdx(blIdx) = blinksIdx(blIdx) - 1;
                        end
                        goOn = 0;
                    end
                end
            else
                goOn = 1;
                while(goOn == 1)
                    if ((samples_x(blinksIdx(blIdx))) - samples_x(blinksIdx(blIdx)+1) >= 2)
                        blinksIdx(blIdx) = blinksIdx(blIdx)+1;
                    else
                        while (samples_x(blinksIdx(blIdx)+ 1) < 0)
                            blinksIdx(blIdx) = blinksIdx(blIdx) + 1;
                        end
                        goOn = 0;
                    end
                end
            end
        end
       
     
     %Delete MS that lie within blinks   
        while (mod(blinks, 2)>0)
            if (trials_work(itrial).left.samples.x(blinksIdx(1)-2) >= 1e07) %single halfblink in beginning
                cut = find(trials_MSextract(itrial).left.Microsaccades.Start < blinksIdx(1));
                
                trials_MSextract(itrial).left.Microsaccades.Start(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.End(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.Merged(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.vPeak(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.DeltaX(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.DeltaY(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.Amplitude(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.Phi(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.StartTime(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.EndTime(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.Duration(cut) = [];
                
                Blink_Indices(1:blinksIdx(1)) = 0; % Write zeros in trials_work at indices that belong to a blink
                
                blinksIdx = blinksIdx(2:blinks);
                blinks = blinks-1;
            else
                cut = find(trials_MSextract(itrial).left.Microsaccades.End > blinksIdx(blinks));
                
                trials_MSextract(itrial).left.Microsaccades.Start(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.End(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.Merged(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.vPeak(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.DeltaX(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.DeltaY(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.Amplitude(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.Phi(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.StartTime(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.EndTime(cut) = [];
                trials_MSextract(itrial).left.Microsaccades.Duration(cut) = [];
                
                Blink_Indices(blinksIdx(end):end) = 0;
                
                blinksIdx = blinksIdx(1:blinks-1);
                blinks = blinks-1;
            end
        end
            
        iiblIdx = 2;
        for ii = 1:blinks/2
            cut1 = find(trials_MSextract(itrial).left.Microsaccades.Start > blinksIdx(iiblIdx-1) & trials_MSextract(itrial).left.Microsaccades.End < blinksIdx(iiblIdx));
            cut2 = find(trials_MSextract(itrial).left.Microsaccades.Start < blinksIdx(iiblIdx-1) & trials_MSextract(itrial).left.Microsaccades.End > blinksIdx(iiblIdx-1));
            cut3 = find(trials_MSextract(itrial).left.Microsaccades.Start < blinksIdx(iiblIdx) & trials_MSextract(itrial).left.Microsaccades.End > blinksIdx(iiblIdx));
            
            cut = [cut1 cut2 cut3];
            trials_MSextract(itrial).left.Microsaccades.Start(cut) = [];
            trials_MSextract(itrial).left.Microsaccades.End(cut) = [];
            trials_MSextract(itrial).left.Microsaccades.Merged(cut) = [];
            trials_MSextract(itrial).left.Microsaccades.vPeak(cut) = [];
            trials_MSextract(itrial).left.Microsaccades.DeltaX(cut) = [];
            trials_MSextract(itrial).left.Microsaccades.DeltaY(cut) = [];
            trials_MSextract(itrial).left.Microsaccades.Amplitude(cut) = [];
            trials_MSextract(itrial).left.Microsaccades.Phi(cut) = [];
            trials_MSextract(itrial).left.Microsaccades.StartTime(cut) = [];
            trials_MSextract(itrial).left.Microsaccades.EndTime(cut) = [];
            trials_MSextract(itrial).left.Microsaccades.Duration(cut) = [];
            
            Blink_Indices(blinksIdx(iiblIdx-1):blinksIdx(iiblIdx)) = 0;
            
            iiblIdx = iiblIdx + 2;
        end
        
        % Replace blink values with interpolation if (blinks >1)
        iiblIdx = 2; %Index of end of each blink
        while (iiblIdx <= blinks)
            tmpIdx = blinksIdx(iiblIdx-1):blinksIdx(iiblIdx);
            trials_MSextract(itrial).left.samples.x(tmpIdx) = linspace(samples_x(blinksIdx(iiblIdx-1)-1), samples_x(blinksIdx(iiblIdx)+1), length(tmpIdx));
            trials_MSextract(itrial).left.samples.y(tmpIdx) = linspace(samples_y(blinksIdx(iiblIdx-1)-1), samples_y(blinksIdx(iiblIdx)+1), length(tmpIdx));
            iiblIdx = iiblIdx +2;
            
        end
    end

    trials_MSextract(itrial).left.samples.Blink_Indices = Blink_Indices;
    
end

clear('samples_x', 'samples_y', 'pupil', 'sortIdx_x', 'sortIdx_y', 'sortIdx_p', 'diffs_sorted', 'blIdx', 'blinksIdx_y', 'blinksIdx_p', 'M', 'Blink_Indices', 'tmpIdx');


%% Calculate measures regarding MS
    
%Number of MS & Amplitude
nMS_pTrial = NaN(1, ntrials_tot);
meanAmp_pTrial = NaN(1, ntrials_tot);
meanXY_pTrial = NaN(2, ntrials_tot);

for itrial = 1:length(trials_MSextract)
    nMS_pTrial(itrial) = size(trials_MSextract(itrial).left.Microsaccades.Start, 2);
    meanAmp_pTrial(itrial) = median(trials_MSextract(itrial).left.Microsaccades.Amplitude);
    meanXY_pTrial(1, itrial) = median(trials_MSextract(itrial).left.Microsaccades.DeltaX);
    meanXY_pTrial(2, itrial) = median(trials_MSextract(itrial).left.Microsaccades.DeltaY);
end


%Number of Saccades
nS_pTrial = NaN(1, ntrials_tot);
% meanAmpS_pTrial = NaN(1, length(trials_work));

for itrial = 1:ntrials_tot
    nS_pTrial(itrial) = size(trials_work(itrial).left.saccade.start, 2);
%     meanAmp_pTrial(itrial) = median(trials_work(itrial).left.saccade.Amplitude);
    
end


%Calculate dispersion ellipse
n_sdH_pTrial = NaN(1, ntrials_tot);
n_sdV_pTrial = NaN(1, ntrials_tot);

for itrial = 1:ntrials_tot
    n_sdH_pTrial(itrial) = sqrt(sum((trials_MSextract(itrial).left.Microsaccades.DeltaX - meanXY_pTrial(1,itrial)).^2)/(ntrials_tot-1));
    n_sdV_pTrial(itrial) = sqrt(sum((trials_MSextract(itrial).left.Microsaccades.DeltaY - meanXY_pTrial(2,itrial)).^2)/(ntrials_tot-1));
end

area_dispEll = NaN(1, ntrials_tot);
p = NaN(1, ntrials_tot);
k = 1;
for itrial = 1:ntrials_tot, p_tmp = min(corrcoef(trials_MSextract(itrial).left.Microsaccades.DeltaX, trials_MSextract(itrial).left.Microsaccades.DeltaY)); p(itrial) = p_tmp(1); end %???
area_dispEll = n_sdV_pTrial .* n_sdH_pTrial .* (k * pi * 2) .* sqrt(1 - p.^2);


%% Calc 2D Velocity Space for moving window of 5 samples
    
deltaT = 1000/250; %?????
thresholdFac = 6;

for itrial = 9:12
    tmpIdx = find(trials_MSextract(itrial).left.samples.Blink_Indices == 1);
    sample_vec = zeros(2, length(trials_MSextract(itrial).left.samples.x(tmpIdx)));
    sample_vec(1, :) = trials_MSextract(itrial).left.samples.x(tmpIdx);
    sample_vec(2, :) = trials_MSextract(itrial).left.samples.y(tmpIdx);
    
    velo_x = zeros(1, length(trials_MSextract(itrial).left.samples.x(tmpIdx)) - 2*2); %vector that will contain velocity data of each trial
    velo_y = zeros(1, length(trials_MSextract(itrial).left.samples.x(tmpIdx)) - 2*2);
    
    vIdx = 1;
    for n = 3:length(trials_MSextract(itrial).left.samples.x(tmpIdx))-2
        velo_x(vIdx) = (sample_vec(1,n+2) + sample_vec(1,n+1) - sample_vec(1,n-1) - sample_vec(1,n-2))/(deltaT/1000 * 6); %factor 6 for n+/-2
        velo_y(vIdx) = (sample_vec(2,n+2) + sample_vec(2,n+1) - sample_vec(2,n-1) - sample_vec(2,n-2))/(deltaT/1000 * 6);
        vIdx = vIdx+1;
    end
    
    figure, plot(velo_x, velo_y);
%     figure, plot(sample_vec(1,:), sample_vec(2,:));
    
    xThreshold = thresholdFac * sqrt((median(velo_x.^2) - median(velo_x)^2));
    yThreshold = thresholdFac * sqrt((median(velo_y.^2) - median(velo_y)^2));
    
    xCenter = 0;
    yCenter = 0;
    
    theta = 0 : 0.01 : 2*pi;
    x = xThreshold * cos(theta) + xCenter;
    y = yThreshold * sin(theta) + yCenter;
    hold all, plot(x, y, '--', 'LineWidth', 3);
    
end
    
clear('tmpIdx');

%% Mean per Trial over Subjects



%% Compare Conditions

condition_names = {'GaussCirclesMasked_Dot', 'CrossedBulleye', 'Bulleye', 'StatCircles'};

Means = table(condition_names', zeros(4,1), zeros(4,1), zeros(4,1), zeros(4,1),  zeros(4,1), 'VariableNames', {'Condition', 'MeanMS', 'MeanAmp', 'MeanDeltaX', 'MeanDeltaY', 'MeanSaccades'}); 

for icond = 1:4
    cond_tmp = find(Exp_block_ID == icond);
%     Means.Condition(icond) = condition_names(icond);
    Means.MeanMS(icond) = mean(nMS_pTrial(cond_tmp));
    Means.MeanAmp(icond) = mean(meanAmp_pTrial(cond_tmp));
    Means.MeanDeltaX(icond) = mean(meanXY_pTrial(1, cond_tmp));
    Means.MeanDeltaY(icond) = mean(meanXY_pTrial(2, cond_tmp));
    Means.MeanSaccades(icond) = mean(nS_pTrial(cond_tmp));
end





%% Plots
figure, plot(nS_pTrial, 'b')
figure, plot(meanAmp_pTrial, 'g')

figure, plot(1:length(trials_MSextract(itrial).left.Microsaccades.Amplitude(:)), trials_MSextract(itrial).left.Microsaccades.Amplitude(:))

figure, plot(number_blinks)

figure, plot(area_dispEll)

%plot x and y values of single trials (set itrial)
figure, plot(1:length(trials_MSextract(itrial).left.samples.x(:)), trials_MSextract(itrial).left.samples.x(:))
hold all, plot(1:length(trials(itrial).left.samples.x(:)), trials(itrial).left.samples.y(:), 'r')


for itrial = 1:5
    figure, plot(trials_work(itrial).left.samples.x, trials_work(itrial).left.samples.y)
    hold all, plot(trials_work(itrial).left.samples.x(trials_MSextract(itrial).left.Microsaccades.Start), trials_work(itrial).left.samples.y(trials_MSextract(itrial).left.Microsaccades.Start), 'r')
end

for itrial = 1:5
%noBlinks
figure, plot(1:length(trials_work(itrial).left.samples.x(:)), trials_work(itrial).left.samples.x(:))
% hold all, plot(1:length(trials_work(itrial).left.samples.x(:)), trials_work(itrial).left.samples.y(:))
end

% %add reference lines to plots
% for ilines = 1:length(condition_list)
%     y = linspace(0, max(, length(trials_new)*10);
%     x = ones(size(y))*iLines*ntrials_pBlock;
%     hold on, plot(x,y,'-')
% end
















