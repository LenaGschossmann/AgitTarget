
function [data] = subjectAnalysis(SbjNumber,MSdataMain)
init_agittarget

% SbjNumber = '7';
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

trials_work = edfstruct;
for itrial = 1:ntrials_tot
    if  ~isstruct(trials_work(itrial).left)
        trials_work(itrial).left = trials_work(itrial).right;
    end
end

%cut-off at begin and end
cutoff_sec = 2; %Cut-off: first and last 2 seconds (equals 250 samples) of each trial
fs = 250;% sampling rate
cutoff = cutoff_sec * fs;

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

for itrial = 1:ntrials_tot
    
    %Blink_Indices contains 0 at every element that has to be disregarded as it lies in a blink
    Blink_Indices = ones(1, length(trials_work(itrial).left.samples.x)); 
    
    samples_x = trials_work(itrial).left.samples.x(:);
    samples_y = trials_work(itrial).left.samples.y(:);
    %     pupil = trials_work(itrial).left.samples.pupil(:);
    
    % check if there are extremely high values in vector which indicate blinks
    if(sum(samples_x(:) >= 1e07)>0) 
        diffs_x = NaN(1, length(samples_x)-1);
        
        diffs_x(:) = abs(diff(samples_x));
        
        [diffs_sorted, sortIdx_x] = sort(diffs_x, 'descend'); % Indices of the sample points with highest differences in x-values
        

        [M, blinks] = max(abs(diff(diffs_sorted))); %blinks/2 = number of blinks
        blinksIdx = sort(sortIdx_x(1:blinks)); %indices of blinks in samples_x in right chronological order
        
        % For uneven number of blinks, check if halfblink in beginning or
        %end
        halfblink = 0; % 0 = no halfblink | 1 = halfblink in beginning | 2 = halfblink in end
        if (mod(blinks,2)>0)
            %halfblink in beginning
            if(trials_work(itrial).left.samples.x(blinksIdx(1)) >= 1e07 && trials_work(itrial).left.samples.x(blinksIdx(1)+2) < 1e06) 
                halfblink = 1;
                blinksStart = (0:2:blinks);
                blinksStart = blinksStart(2:end);
                blinksEnd = (1:2:blinks);
            %halfblink in end
            else 
                halfblink = 2;
                blinksStart = (1:2:blinks);
                blinksEnd = (0:2:blinks);
                blinksEnd = blinksEnd(2:end);
            end
        else %no halfblink
            halfblink = 0;
            blinksStart = (1:2:blinks);
            blinksEnd = (0:2:blinks);
            blinksEnd = blinksEnd(2:end);
        end
                
        %include area around blinks
        
        %threshold for cut-out area around blinks: throw out each next
        %value that is still bigger than the 95% quantile of all x values
        qntl_x = quantile(samples_x, 0.95); 
        
        for blIdx = 1:blinks
            %start of blink
            if (ismember(blIdx, blinksStart))
                while(samples_x(blinksIdx(blIdx)) > qntl_x)
                    if (blinksIdx(blIdx) == 1)%Safety-bar
                       break;
                    else
                        blinksIdx(blIdx) = blinksIdx(blIdx)-1;  % go stepwise back in time and check slope
                    end
                end
            %end of blink    
            else 
                while((samples_x(blinksIdx(blIdx))) > qntl_x)
                    if (blinksIdx(blIdx) == length(blinksIdx))%Safety-bar
                        break;
                    else
                        blinksIdx(blIdx) = blinksIdx(blIdx)+1; % go stepwise forth in time and check slope
                    end 
                end
            end
        end
       
       %set first blinkIndex to 1 if 0
       if(blinksIdx(1) == 0)
           blinksIdx(1) = 1;
       end
       
       trials_MSextract(itrial).left.samples.Blink_Indices = ones(1, length(samples_x));
       
        %% Exclude MS that lie within blinks,start or end in blinks
        % -> exclude MS that lie between the blinkStart&End
        if(halfblink == 1)
            blinksIdx = [1 blinksIdx];
        elseif (halfblink == 2)
            blinksIdx = [blinksIdx length(samples_x)];   
        end
        
        blinksStart = (1:2:length(blinksIdx));
        blinksEnd = (0:2:length(blinksIdx));
        blinksEnd = blinksEnd(2:end);
        
        cut = 0;
       
        for ii = 1:length(blinksStart)   %detect MS in blinks -> code them with 1
            cut = cut + (trials_MSextract(itrial).left.Microsaccades.Start > blinksIdx(blinksStart(ii)) & trials_MSextract(itrial).left.Microsaccades.Start < blinksIdx(blinksEnd(ii)));
            cut = cut + (trials_MSextract(itrial).left.Microsaccades.End > blinksIdx(blinksStart(ii)) & trials_MSextract(itrial).left.Microsaccades.End < blinksIdx(blinksEnd(ii)));  
            % indices of samples that belong to blink = 1
            trials_MSextract(itrial).left.samples.Blink_Indices(blinksIdx(blinksStart(ii)):blinksIdx(blinksEnd(ii))) = 0;
        end
        
        trials_MSextract(itrial).left.samples.Blink_Indices(samples_x < 0) = 0;
        trials_MSextract(itrial).left.samples.Blink_Indices(samples_x > 10^5) = 0;
        
        cut(trials_MSextract(itrial).left.Microsaccades.Amplitude > 10^5) = 0;
        cut(trials_MSextract(itrial).left.Microsaccades.Amplitude < -10^3) = 0;
        
        Fn = fieldnames(trials_MSextract(itrial).left.Microsaccades);
        cut = find(cut == 1 | cut == 2);
        %safe backup of MS data
        MS_backup = trials_MSextract(itrial).left.Microsaccades;
        for ifn = 1:length(Fn)
            trials_MSextract(itrial).left.Microsaccades.(Fn{ifn})(cut) = NaN;
        end
        
            

% %         %%%%%% OLD!!
% %         %Exclude MS that lie within blinks
% %         blinks_backup = blinks;
% %         blinksIdx_backup = blinksIdx;
% %         Fn = fieldnames(trials_MSextract(itrial).left.Microsaccades);
% %         if(halfblink == 1)
% %             cut = find(trials_MSextract(itrial).left.Microsaccades.Start < blinksIdx(1));
% %             for ifn = 1:length(Fn)
% %                     trials_MSextract(itrial).left.Microsaccades.(Fn{ifn})(cut) = [];
% %             end
% %             Blink_Indices(1:blinksIdx(1)) = 0;
% %             blinksEnd = blinksEnd(2:end);
% %         elseif(halfblink == 2 || halfblink == 0)
% %             cut = find(trials_MSextract(itrial).left.Microsaccades.End > blinksIdx(end));
% %             for ifn = 1:length(Fn)
% %                     trials_MSextract(itrial).left.Microsaccades.(Fn{ifn})(cut) = [];
% %             end
% %             Blink_Indices(blinksIdx(end):end) = 0;
% %             blinksStart = blinksStart(2:end);
% %         end
% %         
% %         clear cut;
% %         
% %         for iblinks = 1:length(blinksStart)
% %             cut_tmp = find(trials_MSextract(itrial).left.Microsaccades.End > blinksIdx(blinksStart(iblinks)) & trials_MSextract(itrial).left.Microsaccades.Start < blinksIdx(blinksEnd(iblinks)));
% %             cut = [cut cut_tmp];
% %             for ifn = 1:length(Fn)
% %                     trials_MSextract(itrial).left.Microsaccades.(Fn{ifn})(cut) = [];
% %             end
% %             Blink_Indices(blinksIdx(blinksStart(iblink)):blinksEnd(iblink)) = 0;
% %         end
%         
% %         while (mod(blinks, 2)>0)
% %             if (halfblink == 1) %single halfblink in beginning
% %                 cut = find(trials_MSextract(itrial).left.Microsaccades.Start < blinksIdx(1));
% %                 
% %                 for fn = fieldnames(trials_MSextract(itrial).left.Microsaccades)'
% %                     trials_MSextract(itrial).left.Microsaccades.(fn{1})(cut) = [];
% %                 end
% %                 Blink_Indices(1:blinksIdx(1)) = 0; % Write zeros in trials_work at indices that belong to a blink
% %                 
% %                 blinkLoc = 1; %means beginning
% %                 blinksIdx = blinksIdx(2:blinks);
% %                 blinks = blinks-1;
% %             else   % single halfblink in the end
% %                 cut = find(trials_MSextract(itrial).left.Microsaccades.End > blinksIdx(blinks));
% %                 
% %                 for fn = fieldnames(trials_MSextract(itrial).left.Microsaccades)'
% %                     trials_MSextract(itrial).left.Microsaccades.(fn{1})(cut) = [];
% %                 end
% %                 
% %                 Blink_Indices(blinksIdx(end):end) = 0;
% %                 blinkLoc = 2; %means end
% %                 blinksIdx = blinksIdx(1:blinks-1);
% %                 blinks = blinks-1;
% %             end
% %         end
%         
% %         %Exclude MS that start or end in blinks
% %         blinks = blinks_backup;
% %         blinksIdx = blinksIdx_backup;
% %         if (blinks > 1)
% %             iiblIdx = 2;
% %             for ii = 1:blinks/2
% %                 cut1 = find(trials_MSextract(itrial).left.Microsaccades.Start > blinksIdx(iiblIdx-1) & trials_MSextract(itrial).left.Microsaccades.End < blinksIdx(iiblIdx));
% %                 cut2 = find(trials_MSextract(itrial).left.Microsaccades.Start < blinksIdx(iiblIdx-1) & trials_MSextract(itrial).left.Microsaccades.End > blinksIdx(iiblIdx-1));
% %                 cut3 = find(trials_MSextract(itrial).left.Microsaccades.Start < blinksIdx(iiblIdx)   & trials_MSextract(itrial).left.Microsaccades.End > blinksIdx(iiblIdx));
% %                 
% %                 cut = [cut1 cut2 cut3];
% %                 for fn = fieldnames(trials_MSextract(itrial).left.Microsaccades)'
% %                     trials_MSextract(itrial).left.Microsaccades.(fn{1})(cut) = [];
% %                 end
% %                 %                 trials_MSextract(itrial).left.Microsaccades.Start(cut) = [];
% %                 %                 trials_MSextract(itrial).left.Microsaccades.End(cut) = [];
% %                 %                 trials_MSextract(itrial).left.Microsaccades.Merged(cut) = [];
% %                 %                 trials_MSextract(itrial).left.Microsaccades.vPeak(cut) = [];
% %                 %                 trials_MSextract(itrial).left.Microsaccades.DeltaX(cut) = [];
% %                 %                 trials_MSextract(itrial).left.Microsaccades.DeltaY(cut) = [];
% %                 %                 trials_MSextract(itrial).left.Microsaccades.Amplitude(cut) = [];
% %                 %                 trials_MSextract(itrial).left.Microsaccades.Phi(cut) = [];
% %                 %                 trials_MSextract(itrial).left.Microsaccades.StartTime(cut) = [];
% %                 %                 trials_MSextract(itrial).left.Microsaccades.EndTime(cut) = [];
% %                 %                 trials_MSextract(itrial).left.Microsaccades.Duration(cut) = [];
% %                 
% %                 Blink_Indices(blinksIdx(iiblIdx-1):blinksIdx(iiblIdx)) = 0;
% %                 
% %                 iiblIdx = iiblIdx + 2;
% %             end
% %         end
% %         
        
% %         % Replace blink values with interpolatisortIdx_xon if (blinks >1)
% %         if(mod(blinks, 2)>0) %when there are halfblinks
% %             if (blinkLoc == 1) % second condition necessary if theres just one halfblink
% %                 iiblIdx = 1;
% %                 tmpIdx = 1:blinksIdx(iiblIdx);
% %                 trials_MSextract(itrial).left.samples.x(tmpIdx) = ones(1, length(tmpIdx)) * samples_x(blinksIdx(iiblIdx)+1);
% %                 trials_MSextract(itrial).left.samples.y(tmpIdx) = ones(1, length(tmpIdx)) * samples_y(blinksIdx(iiblIdx)+1);
% %                 blinksIdx = blinksIdx(2:end);
% %             else
% %                 iiblIdx = length(blinksIdx);
% %                 tmpIdx = blinksIdx(iiblIdx):blinksIdx(end);
% %                 trials_MSextract(itrial).left.samples.x(tmpIdx) = ones(1, length(tmpIdx)) * samples_x(blinksIdx(iiblIdx)-1);
% %                 trials_MSextract(itrial).left.samples.y(tmpIdx) = ones(1, length(tmpIdx)) * samples_y(blinksIdx(iiblIdx)-1);
% %                 blinksIdx = blinksIdx(1:end-1);
% %             end
% %         end
% %         iiblIdx = 2; %Index of end of each blink
% %         while (iiblIdx <= blinks)
% %             tmpIdx = blinksIdx(iiblIdx-1):blinksIdx(iiblIdx);
% %             if (blinksIdx(iiblIdx-1) == 1) % case that saccade starts at very beginning -> there is no preceding value
% %                 trials_MSextract(itrial).left.samples.x(tmpIdx) = ones(1, length(tmpIdx)) * samples_x(blinksIdx(iiblIdx)+1);
% %                 trials_MSextract(itrial).left.samples.y(tmpIdx) = ones(1, length(tmpIdx)) * samples_y(blinksIdx(iiblIdx)+1);
% %                 iiblIdx = iiblIdx +2;
% %             elseif (blinksIdx(iiblIdx) == length(blinksIdx)) % case that saccade ends at very end -> there is no following value
% %                 trials_MSextract(itrial).left.samples.x(tmpIdx) = ones(1, length(tmpIdx)) * samples_x(blinksIdx(iiblIdx -1)-1);
% %                 trials_MSextract(itrial).left.samples.y(tmpIdx) = ones(1, length(tmpIdx)) * samples_y(blinksIdx(iiblIdx-1)-1);
% %                 iiblIdx = iiblIdx +2;
% %             else
% %                 trials_MSextract(itrial).left.samples.x(tmpIdx) = linspace(samples_x(blinksIdx(iiblIdx-1)-1), samples_x(blinksIdx(iiblIdx)+1), length(tmpIdx));
% %                 trials_MSextract(itrial).left.samples.y(tmpIdx) = linspace(samples_y(blinksIdx(iiblIdx-1)-1), samples_y(blinksIdx(iiblIdx)+1), length(tmpIdx));
% %                 iiblIdx = iiblIdx +2;
% %             end
% %         end
    end
    
end

clear('samples_x', 'samples_y', 'sortIdx_x', 'diffs_x', 'diffs_sorted', 'blIdx', 'blinksIdx', 'M', 'Blink_Indices', 'tmpIdx');


%% Table with all MS per subject listed

condition_names = {'GaussCirclesMasked_Dot', 'CrossedBulleye', 'Bulleye', 'StatCircles'};

if nargin>1 && MSdataMain
    % benes kram
    data = table();
    for itrial = 1:length(trials_MSextract)
        ms = trials_MSextract(itrial).left.Microsaccades; % all 
        t = struct2table(structfun(@(x)double(x)',ms,'UniformOutput',0)); %transform in double precision values
        t.trial = repmat(itrial,size(t,1),1);
        t.condition = repmat(condition_names(experimentmat.Exp_block_ID(itrial)),size(t,1),1);
        t.subject = repmat(SbjNumber,size(t,1),1);
        data = [data;t];
    end
%     Means  = data; % just to keep your function working...
end


%% Mean MS measures
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

Means = table(condition_names', zeros(4,1), zeros(4,1), zeros(4,1), zeros(4,1),  zeros(4,1), 'VariableNames', {'Condition', 'MeanMS', 'MeanAmp', 'MeanDeltaX', 'MeanDeltaY', 'MeanSaccades'});

for icond = 1:4
    cond_tmp = find(experimentmat.Exp_block_ID== icond);
    %     Means.Condition(icond) = condition_names(icond);
    Means.MeanMS(icond) = mean(nMS_pTrial(cond_tmp));
    Means.MeanAmp(icond) = mean(meanAmp_pTrial(cond_tmp));
    Means.MeanDeltaX(icond) = mean(meanXY_pTrial(1, cond_tmp));
    Means.MeanDeltaY(icond) = mean(meanXY_pTrial(2, cond_tmp));
    Means.MeanSaccades(icond) = mean(nS_pTrial(cond_tmp));
end

if nargin>1 && MSdataMain == 0
    data = Means;
end

end

