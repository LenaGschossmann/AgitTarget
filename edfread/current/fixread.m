function fixmat = fixread( varargin )
% fixmat = fixread( FILENAME, [,'-binocular'] [,'-nocache'] [,'-drift']
%                   [,'-message', MESSAGEFILTER] [,'-subject', SUBJECTINDEX]
%                   [,CALLBACK] )
% 
% WARNING: this version of fixread is for legacy experiment data only
%          use the 'fixations.m' script for new pytrack-style experiment data
error('Legacy fixread version. Please use fixations.m instead. If you must use fixread.m, edit the source file and remove this error manually');
%
% FILENAME
%           FILENAME: String
%           a real filename or a filter mask for filenames
%           Note: FILENAME _MUST_ be an absolute path
% -binocular
%           load data for left and right eye
%           the default is to autoselect depending on calibration results
% -nocache
%           always run edfread on the original EDF
%           the default is to create a .mat file for every .edf file and
%           read from this on repeated invocations. this fails when you
%           change the -message option
% -drift
%           include drift correction results in output
% -message MESSAGEFILTER
%           MESSAGEFILTER: String 
%           look for MESSAGEFILTER messages in EDF to extract image and
%           condition
% -subject SUBJECTINDEX
%           SUBJECTINDEX: Int8 Scalar
%           do not automatically index subjects, but use for all of them
%
% CALLBACK
%           CALLBACK: function pointer
%           method called to convert MESSAGEFILTER messages to [image,
%           condition]
%
% Example:
%
% fixread('filename.edf', 'files/*.edf', '-nocache', '-message', '!V IMGLOAD', @callback)
%
% where IMGLOAD yields eg 'image30condition2.bmp' and you define:
%
% function [image, condition] = callback(msgstr)
%   [mat, ind ] = regexp(msgstr, '\d+', 'match');
%   image = mat(1);
%   condition = mat(2);
% end
% 
% so that [30, 2] = callback('image30condition2.bmp')
%
% TODO : repair the absolute path issue
% 
% $Id: fixread.m 9 2007-06-18 04:46:45Z jsteger $

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % General Layout:
    % initialization of default parameters
    % option parsing and file list creation (filelist())
    % initialization of output struct array (mkeyedata())
    % population of output (addsubject())
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    % image counter for dummy callback
    c_image = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % if no callback function is given this default one is used.
    % MSGSTR contains the file name of the image shown (cf conditionfilter
    % below)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    function [image, condition] = default_callback(msgstr)
        c_image = c_image + 1;
        image   = c_image;
        condition = NaN;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize default parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cb = @default_callback;
    % output both eyes
    binocular = false;
    % output drift correction
    drift = false;
    % read the EDF file from the cache AND write it to cache if it is not cached?
    cache = true;
    % use one subject-index for all files when>0
    fixsubject = 0;
    % the field name which contains the file name of the presented image.
    % this field will be used as a callback function to extract the
    % condition and image index information
    conditionfilter = '!V TRIAL_VAR_DATA';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reads the file names and all the options given during the call
    % of fixread.
    % NOPTS: nargin(=the number of arguments during the call of fixread)
    % RESULT: cell array with the EDF files names we want to include for
    % the conversion. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function result = filelist( nopts )

        % we store the filelist here:    
        result = {};
        function addres(r)
            result{size(result,2)+1} = r;
        end
        f=0;
        while f<nopts % we run over the input arguments!            
            % while loop so we can skip entries on the way for 2 argument
            % parameters
            f = f+1;
            % Fth input argument of fixread
            opt = varargin{f};
            
            % string means filenames or options
            if isstr(opt)
                fdir = dir(opt);
                % got a directory or file:
                if size(fdir, 1) > 0
                    % take the filename or in case if OPT is a file mask
                    % make it run over all the file names
                    % This results in RESULT cell array variable containing
                    % all the filenames
                    for d=(1:size(fdir,1))
                        addres([fileparts(opt), '/', fdir(d).name]);
                    end
                % got options:
                elseif strcmp(opt, '-message')
                    % Messages to be extracted from the EDF.
                    % will be used in image&condition extraction and passed on
                    % to edfread
                    f = f+1;
                    conditionfilter = varargin{f};
                elseif strcmp(opt, '-nocache')
                    cache = false;
                elseif strcmp(opt, '-binocular')
                    binocular = true;
                elseif strcmp(opt, '-drift')
                    drift = true;
                elseif strcmp(opt, '-subject')
                    f = f+1;
                    fixsubject = varargin{f};
                else
                    display(['Ignoring unknown option/non existent file ', 
                            opt, ', parameter ', int2str(f)]);
                end
            % the only function handle we expect is for the callback
            elseif isa(opt, 'function_handle')
                cb = opt;
            else
                display(['Cannot handle argument ',int2str(f),' - ignoring it']);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % parse options, get and check input filenames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    files = filelist(nargin);
    nfiles = size(files, 2);
    if nfiles==0
        display('No valid input files found. Aborting.')
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREFIX will be either 'left', 'right' (-binocular) or '' (default)
    % in the binocular case each eye specific field will contain the
    % first letter of the corresponding eye, e.g. lstart or rstart for
    % the start time of the fixations of the left or right eye
    % respectively.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function mkeyedat(prefix)
        fixmat.([prefix, 'start']) = zeros(0, 0, 'int32');
        fixmat.([prefix, 'end'])   = zeros(0, 0, 'int32');
        fixmat.([prefix, 'x'])     = zeros(0, 0, 'single');
        fixmat.([prefix, 'y'])     = zeros(0, 0, 'single');
        fixmat.([prefix, 'index'])   = zeros(0, 0, 'uint16');
        fixmat.([prefix, 'pupil']) = zeros(0, 0, 'single');
        if drift
            fixmat.([prefix, 'drift']) = zeros(0, 0, 'single');
        end
        % the following will be unique for both eyes once we
        % finished sync'ing of binocular fixations
        fixmat.([prefix, 'subject'])   = zeros(0, 0, 'uint8');
        fixmat.([prefix, 'condition']) = zeros(0, 0, 'uint8');
        fixmat.([prefix, 'image'])     = zeros(0, 0, 'uint16');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize output matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if binocular
        mkeyedat('l');
        mkeyedat('r');
    else
        mkeyedat('');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOW COMES THE CORE PART OF READING THE DATA AND PUTTING INTO THE
    % FIXMAT FORMAT.
    % ADDSUBJECT calls LOADSUBJECT for each edf file and passes the 
    % data on to LOADMONOCULAR or LOADBIMONOCULAR which in turn call
    % LOADBLOCK to add append data to the output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function addsubject(filename, subject)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read the EDF file using edfread OR if CACHE is true and if
        % there is such a cached file (an edf file already transformed
        % to mat) read that one to dont waste time with the conversion
        % If CACHE is true but there is no such a file then it also
        % saves the converted file to the disk. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function [trials, meta] = loadsubject
            cacheedf = [filename,'.mat'];
            
            if cache && size(dir(cacheedf), 1)>0
                load(cacheedf);
                fprintf(1, '.c');
            else
                [trials, meta] = edfread(filename, conditionfilter);
                fprintf(1, '.');
                if cache
                    save(cacheedf, 'trials', 'meta');
                    fprintf(1, '.');
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % adds a chunk of data to the fixmat output
        % eye: a string determining the eye to read data from, 'left'
        %      or 'right'
        % prefix: string. the prefix for eye specific field names
        % sblock: start trial
        % eblock: end trial
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function addblock(eye, prefix, sblock, eblock)
            if sblock<1
                sblock=1;
            end            
            % get the struct containg data for the requested eye and trials
            eyetrials = [trials(sblock:eblock).(eye)];
            
            % run over all trials, adding more info to the struct
            for f=(1:size(eyetrials,2))
                % total number of fixations in this trial
                nfix = size(eyetrials(f).fixation.x, 2);
                % new field: add index of current fixation in trial
                eyetrials(f).fixation.fix = (1:nfix);                               
                % extract the image and condition index from the edf by
                % using the callback
                [image, condition] = cb(trials(f-1+sblock).(conditionfilter).msg);
                % new fields: image, condition
                % store this information for every fixation of this trial
                eyetrials(f).fixation.image = ones(1, nfix)*image;
                eyetrials(f).fixation.condition = ones(1, nfix)*condition;
                % new field: drift
                % if drift flag is on, we also store the the drift
                % correction value done at the end of the trials, again
                % in each fixation
                if drift
                    eyetrials(f).fixation.drift = repmat(eyetrials(f).drift', 1, nfix);
                end
            end
            % EYETRIALS contains many more additional fields now.
            % cat all fixations together
            eyefix = [eyetrials.fixation];
            nfix = size([eyefix.x], 2);
            % append new eyefix to existing output
            fixmat.([prefix,'x'])     = [ [eyefix.x], fixmat.([prefix,'x'])];
            fixmat.([prefix,'y'])     = [ [eyefix.y], fixmat.([prefix,'y'])];
            fixmat.([prefix,'start']) = [ [eyefix.start], fixmat.([prefix,'start'])];
            fixmat.([prefix,'end'])   = [ [eyefix.end], fixmat.([prefix,'end'])];
            fixmat.([prefix,'index'])   = [ [eyefix.index], fixmat.([prefix,'index'])];
            fixmat.([prefix,'pupil']) = [ [eyefix.pupil], fixmat.([prefix,'pupil'])];
            if drift
                fixmat.([prefix,'drift']) = [ [eyefix.drift], fixmat.([prefix,'drift'])];
            end
            fixmat.([prefix,'image'])     = [ eyefix.image, fixmat.([prefix,'image'])];
            fixmat.([prefix,'condition']) = [ eyefix.condition, fixmat.([prefix,'condition'])];
            fixmat.([prefix,'subject'])   = [ ones(1, nfix)*subject, fixmat.([prefix,'subject']) ];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Find the trials where the a calibration is done. Run from the
        % last one toward the first one. In each step add all the trials
        % which are downstream to the current calibration till the end
        % of the experiment or till the next calibration.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function loadmonocular(trials, calib)         
            % number of calibrations realized
            ncalib = size(calib, 1);
            % running from the back, in which trial we calibrated the last time
            lastc  = size(trials, 2)+1;
            % run over calibrations, from the last one to the first one
            for c=(ncalib:-1:1)
                % get current calibration
                cal = calib(c);
                % is the calibration more recent than the last one?
                if cal.trial ~= lastc 
                    % we have a fresh calibration
                    if isstruct(cal.left) && ((cal.left.err_avg <  cal.right.err_avg) || isnan(cal.right.err_avg)) && ( ~ isnan(cal.left.err_avg) )
                        % left eye has better calibration
                        addblock('left', '', cal.trial, lastc-1);
                        lastc = cal.trial;
                    elseif ~ isnan(cal.right.err_avg)
                        % right eye rocks
                        addblock('right', '', cal.trial, lastc-1);
                        lastc = cal.trial;
                    end % no suitable calibration. we will skip this
                end
            end
            if lastc~=0
                display(['Error: no calibration before trial ', int2str(lastc), ' - using left eye']);
                addblock('left', '', 1, lastc);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % just call addblock and add all the information from both eyes
        % todo: sync left and right fixations, adding dummy NaN to the
        % other eye when only one eye fixated
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function loadbinocular(trials)
            ntrials = size(trials, 2);
            addblock('left', 'l', 1, ntrials);
            addblock('right', 'r',1, ntrials);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % addsubject body
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        [trials, meta] = loadsubject();
        if binocular
            loadbinocular(trials)
        else
            loadmonocular(trials, meta.calib)
        end
        fprintf(1, '.');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % main loop over all input files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(1, 'Progress in %%: ');
    for nf=(1:nfiles)
        if fixsubject==0
            addsubject(files{nf}, nf);
            % reset image index for dummy callback
            c_image = 0;
        else
            addsubject(files{nf}, fixsubject);
        end
        fprintf(1, '%1.0d', int8(nf/nfiles*100));
    end
    fprintf(1, ' - done\n');
end
