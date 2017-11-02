
function draw_Cross(screenRes, window, windowRect, targetVisAngle, px_per_deg, duration)

screenFR = Screen('NominalFrameRate', window);
[window_width, window_height] = Screen('WindowSize', window);

% time1 = GetSecs();
% time2 = GetSecs();
% times = [];
% flipTimestamp = clock;

%Set black/white & gray
white = WhiteIndex(0);
black = BlackIndex(0);
colorFac = white/2;
background = 1*colorFac;

%Line parameters

radL = round(0.5*targetVisAngle * px_per_deg); % diameter of outer circle
diamL = 2*radL+1;

radS = round(0.5*(targetVisAngle/3) * px_per_deg);
diamS = 2*radS+1; % diameter of inner circle (degrees)
cross_thickness = 2;

% line1(1:cross_thickness, 1:diamL) = white;
% line2(1:diamL, 1:cross_thickness) = white;
% TexLine1 = Screen('MakeTexture', window, line1);
% TexLine2 = Screen('MakeTexture', window, line2);

bigGrid(1:window_height,1:window_width) = background;
ctr_Grid = size(bigGrid)/2;
% ctr_Grid(2) = round(ctr_Grid(2)/2);
bigGrid(:,:) = background;

% destRect = [0 0 (diamL) (diamL)];
% [destRect, dh, dv] = CenterRect(destRect, windowRect);

%Draw background
TexBG = Screen('MakeTexture', window, bigGrid);
Screen('DrawTexture', window, TexBG);

%Draw Cross
Screen('DrawLine', window, white, ctr_Grid(2)-radL, ctr_Grid(1),ctr_Grid(2)+radL, ctr_Grid(1), cross_thickness);
Screen('DrawLine', window, white, ctr_Grid(2), ctr_Grid(1)-radL,ctr_Grid(2), ctr_Grid(1)+radL, cross_thickness);


Screen('Flip', window);
Screen('Close', TexBG);
WaitSecs(duration);



end