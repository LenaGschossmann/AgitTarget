PsychDefaultSetup(2)


%Settings
%viewing distance [mm]    
dist = 600;   

%% CommandWindowTalk

disp('\n\n SCREEN DISTANCE SHOULD BE 60cm!')
 
%settings
subject_id = input('\n subjectID: ');
holdorfix = input('\n Between-Subject Condition: fixate/hold [f/h]: ', 's');
        if (strcmp(holdorfix, 'f') == true)
            condition_bS = 1;
        elseif(strcmp(holdorfix, 'h')== true)
            condition_bS = 2;
        end
    
%% Set screen & window settings

screen = max(Screen('Screens')); %Experiment-screen
screenCntrl = min(Screen('Screens'));
Screen('Preference', 'SkipSyncTests', 0);
Screen('Preference', 'VisualDebugLevel', 1);
screenRes = Screen('Resolution', screen);

windowCoordinates = [0,0, screenRes.width, screenRes.height];
windowCoordinates = [10,10,2000, 1000];
[window, windowRect] = Screen('openWindow', screen, [],windowCoordinates, [], [], [], 8);   
[window_width, window_height] = Screen('WindowSize', window);
Screen('BlendFunction', window,  GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

screenFR = Screen('NominalFrameRate', window);
 
Screen('ColorRange', window, 255); 
[cx, cy] = RectCenter(windowRect); % centerxx

[dispWidth dispHeight] = Screen('DisplaySize', screen);
px_per_deg = screenRes.height/(2*atan(dispHeight/(2*dist))*(180/pi)); %as the dispHeight is in px and calculated (by PTB as 2.835 pixels/mm) for a subject-monitor distance of 600mm!!


%% Keyboard Settings

[keyboardIndices, productNames, allInfos] = GetKeyboardIndices();
% keyboardIndex=keyboardIndices(strcmp(productNames,'DELL Dell USB Entry Keyboard')); %??
KbName('UnifyKeyNames');  %Enable unified mode of KbName, so KbName accepts identical key names on all operating systems:
spaceKey = KbName('SPACE');


%% Eyelink Eyetracker

resultpath='/home/experiment/experiments/AgitTarget/Experiment/Data';
addpath(genpath(resultpath)) %results path

% ppdev_mex('Close', 1);

%eyetracking
eyetracking=false; %set true
calibrate_eyelink = false; %set true

% if screenFR~=120
%     error('WRONG MONITOR SETTINGS. CHECK FREQUENCY!')
% end

%frameCorrectionTime = 0.006; %in seconds (warum nicht 8ms?)

%setup eyetracker
if eyetracking==1
    EyelinkInit()
    el=EyelinkInitDefaults(window);% win -> PTB window
end

%set subject EDF-filename
while true;
    EDF_name = sprintf('ET_sb_%u.EDF',subject_id);
    %prevent overwriting existing files
    if exist([resultpath EDF_name], 'file')
        ok = input(sprintf('Filename %s already exists, do you want to overwrite it (y/n)?', EDF_name), 's');
        if(strcmp(ok, 'y'))
            break;
        else
            subject_id = input('\n subjectID: ');
        end
    else
        break;
    end
end

%setup eyelink
if eyetracking==1
    setup_eyetracker;
    
    %open log file
    OpenError = Eyelink('OpenFile', EDF_name);
    if OpenError
        error('Eyelink command OpenFile failed (Error: %d)', OpenError)
    end
    
    sessionInfo = sprintf('%s %s','SUBJECTINDEX',num2str(subject_id));
    Eyelink('message','METAEX %s',sessionInfo);
    Eyelink('Message', 'DISPLAY_COORDS %d %d %d %d', 0, 0, window_width, window_height); %***
    
    % disable Matlab keyboard
    ListenChar(2)
    HideCursor();
    
    %start calibration
    if calibrate_eyelink
        fprintf('\n\nEYETRACKING CALIBRATION...');
        EyelinkDoTrackerSetup(el);
        Eyelink('WaitForModeReady', 500); %***
%         [image,map,alpha] = imread(['/home/experiment/experiments/AgitTarget/Experiment/drift_correction.png']); %***
%         fixIndex = Screen('MakeTexture', win.hndl, image); %***
        fprintf('DONE\n\n');
    end
end

if eyetracking==1 && calibrate_eyelink
     Eyelink('message','SYNCTIME'); %this is t=0 for any data following
end

    
%% Experiment Design
 
%condition_bS [1(fixate) 2(hold)]
fixate = {'FIXATE the target\n as still as possible!', 'HOLD your eyes on the target\n as still as possible!'}; 

conditions_wS = [1 2 3 4];

nblocks = 2; %number of blocks per each condition
ntrials_pBlock = 1; %10 
  
  
dur_target = 5;%15; %sec
dur_blank = 3; %sec 
dur_trial = dur_target + dur_blank;

size_target = 0.6; %in degrees
number_circles = 5; %Define number of white&black circles to b e shown


if (length(conditions_wS) > 1)
    
    condition_list = randmat(subject_id);
%     %Pseudorandom design: a b c b c a
%     cond_list_idx = randperm(length(conditions_wS));
%     cond_list_idx = [cond_list_idx fliplr(cond_list_idx)];
%     %swap the middle elements from [a b c c b a] to [a b c b c a]
%     cond_list_idx(length(conditions_wS)+1) = cond_list_idx(length(conditions_wS)+2);
%     cond_list_idx(length(conditions_wS)+2) = cond_list_idx(length(conditions_wS));
%     condition_list = conditions_wS(cond_list_idx); %-> randperm as often as nblocks indicates
else
    condition_list = ones(1,nblocks)*conditions_wS;
end

ntrials_tot = ntrials_pBlock * length(condition_list);

condition_names = {'GaussCirclesMasked_Dot', 'CrossedBulleye', 'Bulleye', 'StatCircles'};
names_list = cell(1, length(condition_list));
for idx = 1:length(condition_list)
    names_list{idx} = condition_names{condition_list(idx)};
end
        
%%Timemanagement
% startTime = 0;         
% stopTime = 0; y

%%Draw blank Background  
white = WhiteIndex(0);
black = BlackIndex(0);
colorFac = white/2; 
background = 1*colorFac; 

blankGrid(1:screenRes.height,1:screenRes.width) = background;
TexBlank = Screen('MakeTexture', window, blankGrid);

%%Text-corner
centerWin = [0 0 screenRes.height/3 screenRes.width/4];
[centerWin, dh, dv] = CenterRect(centerWin, windowRect);
three = '3';
two = '2';  
one = '1'; 

%% Setup variables

Exp_sb_ID = ones(1, ntrials_tot)*subject_id;
Exp_condition_bS = ones(1, ntrials_tot)*condition_bS;

Exp_trial_IDs = 1:ntrials_tot;

Exp_block_trial = repmat(1:ntrials_pBlock, 1, length(condition_list));
Exp_block_ID = reshape(repmat(condition_list, ntrials_pBlock, 1), [1 (ntrials_pBlock*length(condition_list))]);
Exp_block_name = reshape(repmat(names_list, ntrials_pBlock, 1), [1 (ntrials_pBlock*length(condition_list))]);
Exp_blank_times = nan(1, ntrials_tot);

 
%% Experiment blocks

% HideCursor();
trial_id = 1; %unique
trial_pBlock = 1;

%Draw initial information
Screen('TextSize', window, 40);
cond_txt = {'\n\nYour task in the following will be\nto fixate a target as good as you can\nduring the trials.\n\n', '\n\nYour task in the following will be\nto hold your eyes on the target as still as you can\nduring the trials.\n\n'};
welcome = strcat(sprintf('Welcome to our Experiment!\n\n It consists of a total of %d blocks.\n Between the blocks you can take as much time as you need to relax your eyes.', length(condition_list)), cond_txt(condition_bS), 'After each trial you need to wait for three seconds,\nbefore you can start the next trial by pressing Space.');
welcome = char(welcome);

Screen('DrawTexture', window, TexBlank);
DrawFormattedText(window, welcome, 'center', 'center', [255 255 255], [], [], [], 1.5);
Screen('Flip', window);

spacePress(spaceKey);

for block = 1:length(condition_list)
    
    %Draw block information 
    Screen('TextSize', window, 40);
    if (block == 1)
        info = sprintf('This is block number %d.\nPlease press Space when you are ready!', block);
    else
        info = sprintf('This is block number %d.\n Before it starts we will do a short calibration.\n\nPlease press Space when you are ready!', block);
    end
    Screen('DrawTexture', window, TexBlank);
    DrawFormattedText(window, info, 'center', 'center', [255 255 255], [], [], [], 1.5);
    Screen('Flip', window);
 
    spacePress(spaceKey)
    
    %recalibrate after each block
    if eyetracking==1
        if block>1
            if calibrate_eyelink
                fprintf('\n\nEYETRACKING CALIBRATION...')
                EyelinkDoTrackerSetup(el);
                Eyelink('WaitForModeReady', 500);
                fprintf('DONE\n\n')
            end
        end
    end
    
    %Draw Fixating Instruction
    Screen('TextSize', window, 50);
    Screen('DrawTexture', window, TexBlank);
    DrawFormattedText(window, fixate{condition_bS}, 'center', 'center', [255 255 255], [], [], [], 1.5);
    Screen('Flip', window);
    WaitSecs(3);

    
    for trial = 1:ntrials_pBlock  
            
        %start ET recording
        if eyetracking==1 && calibrate_eyelink
            
            %drift correction wanted???
            %drift corection before each trial
%             Screen('DrawTexture', win.hndl, fixIndex);
%             Screen('Flip', win.hndl);
%             Eyelink('WaitForModeReady', 500);
%             EyelinkDoDriftCorrect2(win.el,win.res(1)/2,win.res(2)/2,0)       

            %metadata
            Eyelink('message','TRIALID %d', trial_id); % trial number   
            Eyelink('Command','record_status_message "Block %d/%d, Trial %d/%d"', block, length(condition_list), trial, ntrials_pBlock); % message for the eyetracking recording pc %***
            trialInfo = sprintf('%s %i','trial',trial);
            Eyelink('message','METATR %s', trialInfo);
            trialInfo = sprintf('%s %s','target condition',names_list{block});
            Eyelink('message','METATR %s', trialInfo);
            
            %Start Recording
            Eyelink('StartRecording');
            Eyelink('WaitForModeReady', 100);
        end
        
        %Draw target
        if(condition_list(block) == 1)
            lastFlip = draw_agitCircles_final(screenRes, window, windowRect, size_target, px_per_deg, dur_target, number_circles, 1, 1, 1); %centerDot, gaussCircle, gaussMask
        elseif(condition_list(block) == 2)     
            draw_crossedBulleye(screenRes, window, windowRect, size_target, px_per_deg, dur_target);  %crossed Bulleye
        elseif(condition_list(block) == 3)
            draw_Bulleye(screenRes, window, windowRect, size_target, px_per_deg, dur_target);  %Bulleye
        elseif(condition_list(block) == 4) 
            draw_statCircles(screenRes, window, windowRect, size_target, px_per_deg, dur_target, number_circles, 0, 0); %statCircles
        end
        
        %Unused conditions
%         lastFlip = draw_agitCircles_final(screenRes, window, windowRect, size_target, px_per_deg, dur_target, number_circles, 0, 0, 1);% sharpCircles, gaussMask
%         draw_Cross(screenRes, window, windowRect, size_target, px_per_deg, dur_target); %Cross
%         lastFlip = draw_agitCircles_final(screenRes, window, windowRect, size_target, px_per_deg, dur_target, number_circles, 1, 1, 0); %gaussCircles
        
        %stop ET recording
        if eyetracking==1 && calibrate_eyelink
            Eyelink('StopRecording');
        end
        
        Screen('DrawTexture', window, TexBlank);
        [blankvbl blankStimOnsetTime blankTimestamp] = Screen('Flip', window);

        blank_start = GetSecs();
        
        %Countdown
        Screen('TextSize', window, 150);
        Screen('DrawTexture', window, TexBlank);
        DrawFormattedText(window, three, 'center', 'center', [80 80 80]);
        Screen('Flip', window);
        WaitSecs(1); 
        Screen('DrawTexture', window, TexBlank);
        DrawFormattedText(window, two, 'center', 'center', [80 80 80]);
        Screen('Flip', window);
        WaitSecs(1);
        Screen('DrawTexture', window, TexBlank);
        DrawFormattedText(window, one, 'center', 'center', [80 80 80]);
        Screen('Flip', window); 
        WaitSecs(1);
        Screen('DrawTexture', window, TexBlank);
        Screen('Flip', window);

        if(trial < ntrials_pBlock)
            spacePress(spaceKey)
        end
        
        Exp_blank_times(trial_id) = GetSecs() - blank_start;
        
        trial_id = trial_id +1;
            
%             Screen('FillRect', window, [178 34 34], [0 0 screenRes.width screenRes.height]);
%             Screen('Flip', window);
%             WaitSecs(2);
    end
        
    
end  

disp(condition_list);
disp(Exp_blank_times');

cnames = {'Subject_ID', 'Condition_bS', 'Trial', 'Block_ID', 'Block_Name', 'Trial_perBlock', 'Blank_Times'}
finaltable = table(Exp_sb_ID', Exp_condition_bS', Exp_trial_IDs', Exp_block_ID', Exp_block_name', Exp_block_trial', Exp_blank_times', 'RowNames', cellfun(@int2str, num2cell(Exp_trial_IDs), 'UniformOutput', false), 'VariableNames', cnames);
writetable(finaltable, strcat(resultpath, sprintf('/AgitT_sb_%u.csv', subject_id)));

outputname = [resultpath sprintf('/AgitT_sb_%u', subject_id)];
fullmatfile = sprintf('%s.mat', outputname);

 if eyetracking==1 && calibrate_eyelink
    full_EDF = sprintf('%s.EDF',outputname);
    Eyelink('CloseFile');
    Eyelink('WaitForModeReady', 500);
    Eyelink('ReceiveFile',sprintf('ET_sb_%u.EDF',subject_id),full_EDF);
    Eyelink('WaitForModeReady', 500);
end

save(fullmatfile);
fprintf('Experiment done.\n');


clear screen;
Screen(window, 'Close');
Screen('CloseAll');



 % A known issue: Eyelink('Shutdown') crashing Matlab in 64-bit Linux
    % cf. http://tech.groups.yahoo.com/group/psychtoolbox/message/12732
    if ~IsLinux(true), Eyelink('Shutdown'); end
    % belt_shutdown(belt)
    Screen('CloseAll');
%     Screen('Preference','Verbosity', prevVerbos); % restore previous verbosity
%     Screen('Preference','VisualDebugLevel', prevVisDbg);% restore prev vis dbg
    ListenChar(1) % restore MATLAB keyboard listening (on command window)
    ShowCursor()
    
  
 

       
