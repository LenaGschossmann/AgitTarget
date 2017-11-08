
%initialize with defaults, tell the Eyetracking on what screen to draw
el=EyelinkInitDefaults(window);
Eyelink('command', 'enable_automatic_calibration = YES'); %automatic calibration

screenInfo = sprintf('screen_pixel_coords = 0 0 %d %d',screenRes.width,screenRes.height);
Eyelink('command',screenInfo);

screenInfo = 'screen_phys_coords = -358 202 358 -202'; %physical screen size: Asus
Eyelink('command',screenInfo);

screenInfo= 'marker_phys_coords = -368 212 -368 -212 368 212 368 -212';
Eyelink('command',screenInfo);

% make sure that we get gaze data from the Eyelink and that we have the
% right information in the EDF file
Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS');
Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,BUTTON');
Eyelink('command', 'link_sample_filter = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');

Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS');
Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('command', 'file_sample_filter = LEFT,RIGHT,GAZE,AREA,GAZERES,STATUS');


% set parser (conservative saccade thresholds)
Eyelink('command', 'select_parser_configuration = 0');

%%%%%%%%%% set eyelink calibration %%%%%%%%%%%%%%
% Calibration
% overwrite tracker config

Eyelink('command','generate_default_targets = NO');

Eyelink('command','calibration_colors = 255,255,255 0,0,0');
Eyelink('command','target_size = 1');
el.foregroundcolour=255;
el.backgroundcolour=127;
%el.backgroundcolour=100;
el.msgfontcolour=100;
el.calibrationtargetsize=0.6; %the larget this is, the larger the overall calibration target
el.calibrationtargetwidth=0.11; %the larger this number, the smaller the white dot in the center
el.calibrationtargetcolour=[255 255 255];


EyelinkUpdateDefaults(el);

% Set custom calibration options
% layout on screen for 9 points:
% 5 1 6
% 3 0 4
% 7 2 8

% the pixel coordinates
% left, middle, right
width = window_width;
height = window_height;

l = round( width*0.33);
lm = round(width*0.415);
m = round(width*0.5);
mr = round(width*0.585);
r = round( width*0.66);
% top, center, bottom
t = round( height*0.33);
tc = round(height*0.415);
c = round(height*0.5);
cb = round(height*0.585);
b = round( height*0.66);

Eyelink('command', 'calibration_type = HV5');
% Eyelink('command', 'validation_type = HV5');

% set coordinates for calibration and validation with 5 points
Eyelink('command','calibration_targets = %d,%d %d,%d %d,%d %d,%d %d,%d', m,c, lm,tc, mr,tc, lm,cb, mr,cb);
% % set coordinates for validation
Eyelink('command','validation_targets  = %d,%d %d,%d %d,%d %d,%d %d,%d', m,c, lm,tc, mr,tc, lm,cb, mr,cb);






% 	Calibration Code Bene

% 	#% % Set custom calibration options
% 	#% % layout on screen for 13 points:
% 	#% % 5    1    6
% 	#% %   9    10
% 	#% % 3    0    4
% 	#% %   11   12
% 	#% % 7    2    8

% 	#% left, left-middle, middle, middle right right
% 	scr_w = exp.get('width')
% 	scr_h = exp.get('height')
% 	l = scr_w/2 + (-2*scr_w/8);
% 	lm = scr_w/2 + (-1*scr_w/8);
% 	m = scr_w/2 + (0*scr_w/8);
% 	mr = scr_w/2 + (1*scr_w/8);
% 	r = scr_w/2 + (2*scr_w/8);
% 	#% % top, top-center, center, center-bottom, bottom
% 	t  =  scr_h/2 + (-2*scr_h/8);
% 	tc = scr_h/2 + (-1*scr_h/8);
% 	c  = scr_h/2 + (0*scr_h/8);
% 	cb = scr_h/2 + (1*scr_h/8);
% 	b  = scr_h/2 + (2*scr_h/8);
% 	
% 	# We use only 5 points.
% 	calibTargets = 'calibration_targets  = %d,%d %d,%d %d,%d %d,%d %d,%d'%(m,c, lm,tc, lm,cb, mr,tc, mr,cb);
% 	validTargets = 'validation_targets   = %d,%d %d,%d %d,%d %d,%d %d,%d'%(m,c, lm,tc, lm,cb, mr,tc, mr,cb);







