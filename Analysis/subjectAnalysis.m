
function [Means] = subjectAnalysis(SbjNumber)
    
addpath('/net/home/student/l/lgschossmann/AgitTarget/edfread/build/linux64');

% SbjNumber = '4';

curSub = sprintf('AgitT_sb_%s', SbjNumber);

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
%     pupil = trials_work(itrial).left.samples.pupil(:);
    
    if(sum(samples_x(:) >= 1e07)>0) %check if there are extremely high values in vector which indicate blinks
        diffs_x = NaN(1, length(samples_x)-1);

        diffs_x(:) = abs(diff(samples_x));

        [diffs_sorted, sortIdx_x] = sort(diffs_x, 'descend'); % Indices of the sample points with highest differences in x-values

        %number of blinks
        [M, blinks] = max(abs(diff(diffs_sorted))); %blinks/2 = number of blinks
        blinksIdx = sort(sortIdx_x(1:blinks)); %indices of blinks in samples_x in right chronological order
        
        %include area around blinks
        for blIdx = 1:blinks
            if (mod(blIdx,2)>0)
                goOn = 1;
                while(goOn == 1)
                    if (blinksIdx(blIdx) == 1)
                        goOn = 0;
                    else
                        if ((samples_x(blinksIdx(blIdx)) - samples_x(blinksIdx(blIdx)-1)) >= 2) % 2 is randomly chosen
                            blinksIdx(blIdx) = blinksIdx(blIdx)-1;
                        else
                            while (samples_x(blinksIdx(blIdx)-1) < 0)
                                blinksIdx(blIdx) = blinksIdx(blIdx) - 1;
                            end
                            goOn = 0;
                        end
                    end
                end
            else
                goOn = 1;
                while(goOn == 1)
                    if (blinksIdx(blIdx) == length(blinksIdx))
                        goOn = 0;
                    else
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
        end
       
     
     %Exclude MS that lie within blinks
     blinks_backup = blinks;
     blinksIdx_backup = blinksIdx;
        while (mod(blinks, 2)>0)
            if (trials_work(itrial).left.samples.x(blinksIdx(1)) >= 1e07 && trials_work(itrial).left.samples.x(blinksIdx(1)+2) < 1e06) %single halfblink in beginning
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
                   
                blinkLoc = 1; %means beginning
                blinksIdx = blinksIdx(2:blinks);
                blinks = blinks-1;         
            else   % single halfblink in the end
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
                blinkLoc = 2; %means end
                blinksIdx = blinksIdx(1:blinks-1);
                blinks = blinks-1;
            end
        end
        
        %Exclude MS that start or end in blinks
        blinks = blinks_backup;
        blinksIdx = blinksIdx_backup;
        if (blinks > 1)
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
        end

        
        % Replace blink values with interpolatisortIdx_xon if (blinks >1)
         if(mod(blinks, 2)>0) %when there are halfblinks
            if (blinkLoc == 1) % second condition necessary if theres just one halfblink
                iiblIdx = 1;
                tmpIdx = 1:blinksIdx(iiblIdx);
                trials_MSextract(itrial).left.samples.x(tmpIdx) = ones(1, length(tmpIdx)) * samples_x(blinksIdx(iiblIdx)+1);
                trials_MSextract(itrial).left.samples.y(tmpIdx) = ones(1, length(tmpIdx)) * samples_y(blinksIdx(iiblIdx)+1);
                blinksIdx = blinksIdx(2:end);
            else
                iiblIdx = length(blinksIdx);
                tmpIdx = blinksIdx(iiblIdx):blinksIdx(end);
                trials_MSextract(itrial).left.samples.x(tmpIdx) = ones(1, length(tmpIdx)) * samples_x(blinksIdx(iiblIdx)-1);
                trials_MSextract(itrial).left.samples.y(tmpIdx) = ones(1, length(tmpIdx)) * samples_y(blinksIdx(iiblIdx)-1);
                blinksIdx = blinksIdx(1:end-1);
            end
         end
         iiblIdx = 2; %Index of end of each blink
         while (iiblIdx <= blinks)
             tmpIdx = blinksIdx(iiblIdx-1):blinksIdx(iiblIdx);
             if (blinksIdx(iiblIdx-1) == 1) % case that saccade starts at very beginning -> there is no preceding value
                 trials_MSextract(itrial).left.samples.x(tmpIdx) = ones(1, length(tmpIdx)) * samples_x(blinksIdx(iiblIdx)+1);
                 trials_MSextract(itrial).left.samples.y(tmpIdx) = ones(1, length(tmpIdx)) * samples_y(blinksIdx(iiblIdx)+1);
                 iiblIdx = iiblIdx +2;
             elseif (blinksIdx(iiblIdx) == length(blinksIdx)) % case that saccade ends at very end -> there is no following value
                 trials_MSextract(itrial).left.samples.x(tmpIdx) = ones(1, length(tmpIdx)) * samples_x(blinksIdx(iiblIdx -1)-1);
                 trials_MSextract(itrial).left.samples.y(tmpIdx) = ones(1, length(tmpIdx)) * samples_y(blinksIdx(iiblIdx-1)-1);
                 iiblIdx = iiblIdx +2;
             else
                 trials_MSextract(itrial).left.samples.x(tmpIdx) = linspace(samples_x(blinksIdx(iiblIdx-1)-1), samples_x(blinksIdx(iiblIdx)+1), length(tmpIdx));
                 trials_MSextract(itrial).left.samples.y(tmpIdx) = linspace(samples_y(blinksIdx(iiblIdx-1)-1), samples_y(blinksIdx(iiblIdx)+1), length(tmpIdx));
                 iiblIdx = iiblIdx +2;
             end
         end
    end

    trials_MSextract(itrial).left.samples.Blink_Indices = Blink_Indices;
    
end

clear('samples_x', 'samples_y', 'sortIdx_x', 'diffs_x', 'diffs_sorted', 'blIdx', 'blinksIdx', 'M', 'Blink_Indices', 'tmpIdx');


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


end


