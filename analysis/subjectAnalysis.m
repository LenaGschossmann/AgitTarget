
function [experimentmat MS_data SAC_data trials_OV cond_OV] = subjectAnalysis(SbjNumber,MSdata, SACdata, overview, cleaned)
init_agittarget

% SbjNumber = '5'; 
% cleaned = 1; %Cleaned = Data without trials with samples too far from Fixpoint
    
if isnumeric(SbjNumber)
    SbjNumber = num2str(SbjNumber);
end

curSub = sprintf('AgitT_sb_%s', SbjNumber);
 
edfpath = fullfile(datapath,'Data',sprintf('%s.EDF',curSub));
matpath = fullfile(datapath,'Data',sprintf('%s.mat',curSub));

experimentmat = load(matpath); % .mat file

[edfstruct, info] = edfread(edfpath, 'write'); % EDF data

ntrials_tot = experimentmat.ntrials_tot;


%% Create working trials structure

trials_raw = edfstruct;
trEye = repmat('l', ntrials_tot, 1);
for itrial = 1:ntrials_tot
    if  ~isstruct(trials_raw(itrial).left)
        trials_raw(itrial).left = trials_raw(itrial).right;
        trEye(itrial) = 'r';
    end
end

%cutMS-off at begin and end
cutoff_sec = 2; %cutMS-off: first and last 2 seconds (equals 250 samples) of each trial
fs = 250;% sampling rate
cutoff = cutoff_sec * fs; 

for itrial = 1:ntrials_tot
    borders = [(cutoff+1) length(trials_raw(itrial).left.samples.x)-cutoff];
    trials_raw(itrial).left.samples.time = trials_raw(itrial).left.samples.time(borders(1):borders(2));
    trials_raw(itrial).left.samples.x = trials_raw(itrial).left.samples.x(borders(1):borders(2));
    trials_raw(itrial).left.samples.y = trials_raw(itrial).left.samples.y(borders(1):borders(2));
    trials_raw(itrial).left.samples.pupil = trials_raw(itrial).left.samples.pupil(borders(1):borders(2));
end


%%...as not all are sampled at same frequency
if SbjNumber == 7
    SampleRate = 500;
else
    SampleRate = 250;
end

%% Extract Microsaccades
[trials_MSextract]= edfExtractMicrosaccades(trials_raw, SampleRate);
% MS parameters
%     --> Amplitude: in px
%     --> Phi: Amplitude in vis deg
%     --> normgx & normgy: distance in px
%     --> vPeak: Peak velocity


%% Disregard all MS during blinks if blinks were not cutMS out before

valid_trials = ones(ntrials_tot, 1);
norm_numberMS = zeros(ntrials_tot,1);
MS_StartX = [];
MS_StartY = [];
MS_EndX = [];
MS_EndY = [];

for itrial = 1:ntrials_tot
    
    %Good_Values contains 0 at every element that has to be disregarded as it lies in a blink
%     Good_Values = ones(1, length(trials_raw(itrial).left.samples.x)); 
    
    samples_x = trials_raw(itrial).left.samples.x(:);
    samples_y = trials_raw(itrial).left.samples.y(:);
    %     pupil = trials_raw(itrial).left.samples.pupil(:);
    
    % here are all values that need to be eliminated in MS and SAC coded
    % with 1
    cutMS = zeros(1, length(trials_MSextract(itrial).left.Microsaccades.Start));
    cutSAC = zeros(1, length(trials_MSextract(itrial).left.saccade.start)); %saccades
    %and here are all values that are okay in samples coded with 1
    trials_MSextract(itrial).left.samples.Good_Values = ones(1, length(samples_x));
    
    % check if there are extremely high values in vector which indicate blinks
    if(sum(samples_x(:) >= 1e07)>0) 
        diffs_x = NaN(1, length(samples_x)-1);
        
        diffs_x(:) = abs(diff(samples_x));
        
        [diffs_sorted, sortIdx_x] = sort(diffs_x, 'descend'); % Indices of the sample points with highest differences in x-values

        [M, blinks] = max(abs(diff(diffs_sorted))); %blinks/2 = number of blinks
        blinksIdx = sort(sortIdx_x(1:blinks)); %indices of blinks in samples_x in right chronological order
        
                
%         %%include area around blinks
%         %threshold for cutMS-out area around blinks: throw out each next
%         %value that is still bigger than the 95% quantile of all x values
%         qntl_x = quantile(samples_x, 0.95); 
%         
%         for blIdx = 1:blinks
%             %start of blink
%             if (ismember(blIdx, blinksStart))
%                 while(samples_x(blinksIdx(blIdx)) > qntl_x)
%                     if (blinksIdx(blIdx) == 1)%Safety-bar
%                        break;
%                     else
%                         blinksIdx(blIdx) = blinksIdx(blIdx)-1;  % go stepwise back in time and check slope
%                     end
%                 end
%             %end of blink    
%             else 
%                 while((samples_x(blinksIdx(blIdx))) > qntl_x)
%                     if (blinksIdx(blIdx) == length(blinksIdx))%figure, plot(1:length(trials_MSextract(itrial).left.samples.x(:)), trials_MSextract(itrial).left.samples.x(:))Safety-bar
%                         break;
%                     else
%                         blinksIdx(blIdx) = blinksIdx(blIdx)+1; % go stepwise forth in time and check slope
%                     end 
%                 end
%             end
%         end
%        
        
         %% For uneven number of blinks, check if halfblink in beginning or
        %end
        halfblink = 0; % 0 = no halfblink | 1 = halfblink in beginning | 2 = halfblink in end
        if (mod(blinks,2)>0)
            %halfblink in beginning
            if(trials_raw(itrial).left.samples.x(blinksIdx(1)) >= 1e07 && trials_raw(itrial).left.samples.x(blinksIdx(1)+2) < 1e06) 
                halfblink = 1;
                blinksIdx = [1 blinksIdx];
            %halfblink in end
            else 
                halfblink = 2;
                blinksIdx = [blinksIdx length(samples_x)];
            end
        else %no halfblink
            halfblink = 0;
        end
        
        blinksStart = (1:2:length(blinksIdx));
        blinksEnd = (0:2:length(blinksIdx));
        blinksEnd = blinksEnd(2:end);
        
       %set first blinkIndex to 1 if 0
       if(blinksIdx(1) == 0)    
           blinksIdx(1) = 1;
       end    


        %% cutMS out 200ms before and 500ms after blink
        blinksIdx(blinksStart) = blinksIdx(blinksStart)-200;
        blinksIdx(blinksEnd) = blinksIdx(blinksEnd)+500;
        %In case of blinkIndex being too close to start or end
        blinksIdx(blinksIdx<0) = 1;
        blinksIdx(blinksIdx>length(samples_x)) = length(samples_x);
        %In case of blinks that are overlapping now
        overlap = find(diff(blinksIdx(2:end)) <= 0)+1;
        blinksIdx(overlap) = NaN;  %blink starts
        blinksIdx(overlap+1) = NaN;

        blinksIdx = blinksIdx(~(isnan(blinksIdx))); %and delete NaNs
        
        if (find(blinksIdx == 1) > 1) %in case you have several elements being 1 now
            blinksIdx = blinksIdx(~(find(blinksIdx == 1)));
            blinksIdx = [1 blinksIdx];
        end
        if (sum(blinksIdx == length(samples_x)) > 1) %...or several elements being the last value
            blinksIdx(blinksIdx == length(samples_x)) = [];
            blinksIdx = [blinksIdx length(samples_x)];
        end
                    
       
        %% Exclude MS and SAC that lie within blinks,start or end in blinks
        % -> exclude MS that lie between the blinkStart&End

        blinksStart = (1:2:length(blinksIdx));
        blinksEnd = (0:2:length(blinksIdx));
        blinksEnd = blinksEnd(2:end);
       
        for ii = 1:length(blinksStart)   %detect MS in blinks -> code them with 1
            cutMS = cutMS + (trials_MSextract(itrial).left.Microsaccades.Start >= blinksIdx(blinksStart(ii)) & trials_MSextract(itrial).left.Microsaccades.Start <= blinksIdx(blinksEnd(ii)));
            cutMS = cutMS + (trials_MSextract(itrial).left.Microsaccades.End >= blinksIdx(blinksStart(ii)) & trials_MSextract(itrial).left.Microsaccades.End <= blinksIdx(blinksEnd(ii)));  
            %and saccades
            cutSAC = cutSAC + (trials_MSextract(itrial).left.saccade.start >= blinksIdx(blinksStart(ii)) & trials_MSextract(itrial).left.saccade.start <= blinksIdx(blinksEnd(ii)));
            cutSAC = cutSAC + (trials_MSextract(itrial).left.saccade.end >= blinksIdx(blinksStart(ii)) & trials_MSextract(itrial).left.saccade.end <= blinksIdx(blinksEnd(ii))); 
            % indices of samples that belong to blink = 0
            trials_MSextract(itrial).left.samples.Good_Values(blinksIdx(blinksStart(ii)):blinksIdx(blinksEnd(ii))) = 0;
        end
    end
        
        %% Exclude impossible values --> cut out MS that contain samples that are marked with 0 in Good_Values
        badVals = trials_MSextract(itrial).left.samples.x > experimentmat.screenRes.width | trials_MSextract(itrial).left.samples.x < 0;
        badVals = badVals + (trials_MSextract(itrial).left.samples.y > experimentmat.screenRes.height | trials_MSextract(itrial).left.samples.y < 0);
        trials_MSextract(itrial).left.samples.Good_Values(badVals > 0) = 0; %mark impossible values in Good_Values
        %check if there are badValues lying within a blink
        for ii = 1:length(trials_MSextract(itrial).left.Microsaccades.Start)
            if (sum(find(badVals > 0) >= trials_MSextract(itrial).left.Microsaccades.Start(ii) & find(badVals > 0) <= trials_MSextract(itrial).left.Microsaccades.End(ii)) > 0)
                cutMS(ii) = 1;
            end
        end
        
        %Change Microsaccade fields
        Fn = fieldnames(trials_MSextract(itrial).left.Microsaccades);
        cutMS = find(cutMS > 0);
        %safe backup of MS data
        MS_backup = trials_MSextract(itrial).left.Microsaccades;
        for ifn = 1:length(Fn)
            trials_MSextract(itrial).left.Microsaccades.(Fn{ifn})(cutMS) = [];
        end
        
        %do the same for saccades
        for ii = 1:length(trials_MSextract(itrial).left.saccade.start)
            if (sum(find(badVals > 0) >= trials_MSextract(itrial).left.saccade.start(ii) & find(badVals > 0) <= trials_MSextract(itrial).left.saccade.end(ii)) > 0)
                cutSAC(ii) = 1;
            end
        end
        Fn = fieldnames(trials_MSextract(itrial).left.saccade);
        Fn_new = {'Start', 'Sx', 'Sy', 'End', 'Ex', 'Ey', 'Speed'};
        cutSAC = find(cutSAC > 0);
        %safe backup of MS data
        SAC_backup = trials_MSextract(itrial).left.saccade;
        for ifn = 1:length(Fn)
            trials_MSextract(itrial).left.saccade.(Fn{ifn})(cutSAC) = [];
            [trials_MSextract(itrial).left.saccade.(Fn_new{ifn})] = trials_MSextract(itrial).left.saccade.(Fn{ifn});
            trials_MSextract(itrial).left.saccade = rmfield(trials_MSextract(itrial).left.saccade, (Fn{ifn}));
        end        
        
        %Check if any MS amplitude exceeds 2degree
        if ((1/experimentmat.px_per_deg).* trials_MSextract(itrial).left.Microsaccades.Amplitude > 2)
            warning(sprintf('There are Microsaccades in subject %d, trial %d that exceed 2 degrees in their amplitude', SbjNumber, itrial))
        end
        
        % Calculate mean distance of samples per trial from fixation point
        x_dev = samples_x(trials_MSextract(itrial).left.samples.Good_Values == 1)-(experimentmat.screenRes.width/2);
        y_dev = samples_y(trials_MSextract(itrial).left.samples.Good_Values == 1)-(experimentmat.screenRes.height/2);
        mean_amp_xy = mean(sqrt(x_dev.^2 + y_dev.^2));
        if (mean_amp_xy > 5 * experimentmat.px_per_deg) %threshold of x times visual degreen -> when samples lie outside, trial dismissed
            valid_trials(itrial) = 0;
        end
        
        %Dismiss trials with too many blinks (when usable values are less
        %than 25%
        if sum(trials_MSextract(itrial).left.samples.Good_Values) < 0.25*length(trials_MSextract(itrial).left.samples.Good_Values)
            valid_trials(itrial) = 0;
        end
        
        %Normalize Number of MS by valid sample points
        norm_numberMS(itrial) = length(trials_MSextract(itrial).left.Microsaccades.Start) / sum(trials_MSextract(itrial).left.samples.Good_Values);
 
end

clear('samples_x', 'samples_y', 'sortIdx_x', 'diffs_x', 'diffs_sorted', 'blIdx', 'blinksIdx', 'M', 'Good_Values', 'tmpIdx');



%% Mean MS measures
%Number of MS & Amplitude
nMS_pTrial = NaN(1, ntrials_tot);
meanAmp_pTrial = NaN(1, ntrials_tot);
meanXY_pTrial = NaN(2, ntrials_tot);

for itrial = 1:length(trials_MSextract)
    nMS_pTrial(itrial) = size(trials_MSextract(itrial).left.Microsaccades.Start, 2);
    meanAmp_pTrial(itrial) = median(trials_MSextract(itrial).left.Microsaccades.Amplitude);
end


%% Tables

condition_names = {'GaussCirclesMasked_Dot', 'CrossedBulleye', 'Bulleye', 'StatCircles'};

    %% Table with overview over trials and conditions
    if nargin>1 && overview
        trials_OV = table();
        trials_OV.trial =  [1:ntrials_tot]';
        trials_OV.condition = condition_names(experimentmat.Exp_block_ID)';
        trials_OV.trackedEye = trEye;
        trials_OV.subject = repmat(SbjNumber,ntrials_tot,1);
        trials_OV.valid_trial = valid_trials;
        trials_OV.numberMS = nMS_pTrial';
        trials_OV.norm_numberMS = norm_numberMS;
        trials_OV.meanAmpMS = meanAmp_pTrial';

        cond_OV = table();
        cond_OV.condition = condition_names';
        cond_OV.subject = repmat(SbjNumber, length(condition_names), 1);

        for icond = 1:length(condition_names)
            % number of valid trials per condition
            cond_OV.valid_trials(icond) = sum(valid_trials(strcmp(condition_names{icond}, trials_OV.condition) == 1));
            tmpT = strcmp(condition_names{icond}, trials_OV.condition) == 1 & trials_OV.valid_trial == 1; 
            % mean number of MS per condition
            cond_OV.mean_numberMS(icond) = mean(trials_OV.numberMS(tmpT));
            %mean number of normalized number of MS per condition
            cond_OV.mean_norm_numberMS(icond) = mean(trials_OV.norm_numberMS(tmpT));
            %mean MS amplitude per condition
            cond_OV.mean_AmpMS(icond) = mean(trials_OV.meanAmpMS(tmpT));
        end
        clear('tmpT');
    end
    
    
    %% Table with all MS per subject listed - cleaned
    
    if cleaned
        cleanTrials = (find(valid_trials == 1));
        %Microsaccades detected by E&K
        if nargin>1 && MSdata
            MS_data = table();
            for itrial = cleanTrials'
                ms = trials_MSextract(itrial).left.Microsaccades; % all 
                tMS = struct2table(structfun(@(x)double(x)',ms,'UniformOutput',0)); %transform in double precision values
                tMS.subject = repmat(SbjNumber,size(tMS,1),1);
                tMS.trial = repmat(itrial,size(tMS,1),1);
                tMS.condition = repmat(condition_names(experimentmat.Exp_block_ID(itrial)),size(tMS,1),1);
                tMS.trackedEye = repmat(trEye(itrial), size(tMS,1),1);
                tMS.StartX = trials_MSextract(itrial).left.samples.x(trials_MSextract(itrial).left.Microsaccades.Start)';
                tMS.StartY = trials_MSextract(itrial).left.samples.y(trials_MSextract(itrial).left.Microsaccades.Start)';
                tMS.EndX = trials_MSextract(itrial).left.samples.x(trials_MSextract(itrial).left.Microsaccades.End)';
                tMS.EndY = trials_MSextract(itrial).left.samples.y(trials_MSextract(itrial).left.Microsaccades.End)';
                MS_data = [MS_data;tMS];
            end
            MS_data.Start = [];
            MS_data.End = [];
        end

        %Saccades detected by Eyelink
        if nargin>1 && SACdata
            SAC_data = table();
            for itrial = cleanTrials'
                sac = trials_MSextract(itrial).left.saccade; % all 
                tSAC = struct2table(structfun(@(x)double(x)',sac,'UniformOutput',0)); %transform in double precision values
                tSAC.subject = repmat(SbjNumber,size(tSAC,1),1);
                tSAC.trial = repmat(itrial,size(tSAC,1),1);
                tSAC.condition = repmat(condition_names(experimentmat.Exp_block_ID(itrial)),size(tSAC,1),1);
                tSAC.trackedEye = repmat(trEye(itrial), size(tSAC,1),1);
                tSAC.StartX = trials_MSextract(itrial).left.saccade.Sx';
                tSAC.StartY = trials_MSextract(itrial).left.saccade.Sy';
                tSAC.EndX = trials_MSextract(itrial).left.saccade.Ex';
                tSAC.EndY = trials_MSextract(itrial).left.saccade.Ex';
                SAC_data = [SAC_data;tSAC];
            end  
        end

%         tTrial = table();
%         tTrial.subject = repmat(SbjNumber, length(cleanTrials), 1);
%         tTrial.trial = cleanTrials;
%         tTrial.condition = experimentmat.Exp_block_name';
% 
%         for itrial = cleanTrials'
%             tTrial.tnumMS(itrial) = length(trials_MSextract(itrial).left.Microsaccades.Start);
%         end
%         t.trial.trackedEye = trEye;
%         tTrial.Properties.VariableNames = {'subject', 'trial', 'condition', 'numMS'};
    
    
     %% Table with all MS per subject listed - NOT cleaned
     else
        %Microsaccades detected by E&K
        if nargin>1 && MSdata
            MS_data = table();
            for itrial = 1:length(trials_MSextract)
                ms = trials_MSextract(itrial).left.Microsaccades; % all 
                tMS = struct2table(structfun(@(x)double(x)',ms,'UniformOutput',0)); %transform in double precision values
                tMS.trial = repmat(itrial,size(tMS,1),1);
                tMS.condition = repmat(condition_names(experimentmat.Exp_block_ID(itrial)),size(tMS,1),1);
                tMS.trackedEye = repmat(trEye(itrial), size(tMS,1),1);
                tMS.subject = repmat(SbjNumber,size(tMS,1),1);
                MS_data = [MS_data;tMS];
            end
        %     Means  = MS_data; % just to keep your function working...
        end

        %Saccades detected by Eyelink
        if nargin>1 && SACdata
            SAC_data = table();
            for itrial = 1:length(trials_MSextract)
                sac = trials_MSextract(itrial).left.saccade; % all 
                tSAC = struct2table(structfun(@(x)double(x)',sac,'UniformOutput',0)); %transform in double precision values
                tSAC.trial = repmat(itrial,size(tSAC,1),1);
                tSAC.condition = repmat(condition_names(experimentmat.Exp_block_ID(itrial)),size(tSAC,1),1);
                tSAC.trackedEye = repmat(trEye(itrial), size(tSAC,1),1);
                tSAC.subject = repmat(SbjNumber,size(tSAC,1),1);
                SAC_data = [SAC_data;tSAC];
            end
        end


        tTrial = table();
        tTrial.subject = repmat(SbjNumber, length(trials_MSextract), 1);
        tTrial.trial = [1:80]';
        tTrial.condition = experimentmat.Exp_block_name';

        for itrial = 1:length(trials_MSextract)
            tTrial.tnumMS(itrial) = length(trials_MSextract(itrial).left.Microsaccades.Start);
        end

        t.trial.trackedEye = trEye;
        tTrial.Properties.VariableNames = {'subject', 'trial', 'condition', 'numMS'};
     end

       
% 
% %Number of Saccades
% nS_pTrial = NaN(1, ntrials_tot);
% % meanAmpS_pTrial = NaN(1, length(trials_raw));
% 
% for itrial = 1:ntrials_tot
%     nS_pTrial(itrial) = size(trials_raw(itrial).left.saccade.start, 2);
%     %     meanAmp_pTrial(itrial) = median(trials_raw(itrial).left.saccade.Amplitude);
% end
% 
% Means = table(condition_names', zeros(4,1), zeros(4,1), zeros(4,1), zeros(4,1),  zeros(4,1), 'VariableNames', {'Condition', 'MeanMS', 'MeanAmp', 'MeanDeltaX', 'MeanDeltaY', 'MeanSaccades'});
% 
% for icond = 1:4
%     cond_tmp = find(experimentmat.Exp_block_ID== icond);
%     %     Means.Condition(icond) = condition_names(icond);
%     Means.MeanMS(icond) = mean(nMS_pTrial(cond_tmp));
%     Means.MeanAmp(icond) = mean(meanAmp_pTrial(cond_tmp));
%     Means.MeanDeltaX(icond) = mean(meanXY_pTrial(1, cond_tmp));
%     Means.MeanDeltaY(icond) = mean(meanXY_pTrial(2, cond_tmp));
%     Means.MeanSaccades(icond) = mean(nS_pTrial(cond_tmp));
% end

if nargin>1 && MSdata == 0
    MS_data = Means;
end

end












