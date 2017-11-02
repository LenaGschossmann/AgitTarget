
function draw_Bulleye(screenRes, window, windowRect, targetVisAngle, px_per_deg, duration)

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

%Circle & Line parameters

radL = round(0.5*targetVisAngle * px_per_deg); % diameter of outer circle
diamL = 2*radL+1;
radS = round(0.5*0.2 * px_per_deg); %Inner circle with diameter of 0.2°
diamS = 2*radS+1; % diameter of inner circle (degrees)

% line1(1:diamS, 1:diamL) = white;
% line2(1:diamL, 1:diamS) = white;
% TexLine1 = Screen('MakeTexture', window, line1);
% TexLine2 = Screen('MakeTexture', window, line2);

bigGrid(1:window_height,1:window_width) = background;
ctr_Grid = size(bigGrid)/2;
bigGrid(:,:) = background;

%Circles
circleRect = ones(diamL, diamL)*background; %Rectangle to plot circles in

destRect = [0 0 (diamL) (diamL)];
[destRect, dh, dv] = CenterRect(destRect, windowRect);

%Define which fields are lying within the circle and make circle texture
X = repmat((1:diamL), diamL, 1);
Y = repmat((1:diamL)', 1, diamL);
ctr_XY = (diamL+1)/2;

circleRect(sqrt((X-ctr_XY).^2 + (Y-ctr_XY).^2) <= radL) = black;
circleRect(sqrt((X-ctr_XY).^2 + (Y-ctr_XY).^2) <= radS) = white;
TexCircle = Screen('MakeTexture', window, circleRect);

%Draw background
TexBG = Screen('MakeTexture', window, bigGrid);
Screen('DrawTexture', window, TexBG);


%Draw Bulleye
Screen('drawTexture', window, TexCircle, [], destRect);
% Screen('FillOval', window, colorOval, [ctr_Grid(1)-radL, ctr_Grid(2)-radL, ctr_Grid(1)+radL, ctr_Grid(2)+radL], diamL);
% Screen('DrawLine', window, colorCross, ctr_Grid(2)-radL, ctr_Grid(1),ctr_Grid(2)+radL, ctr_Grid(1), diamS);
% Screen('DrawLine', window, colorCross, ctr_Grid(2), ctr_Grid(1)-radL,ctr_Grid(2), ctr_Grid(1)+radL, diamS);
% Screen('FillOval', window, colorCross, [ctr_Grid(1)-radS, ctr_Grid(2)-radS, ctr_Grid(1)+radS, ctr_Grid(2)+radS], diamS);

Screen('Flip', window);
Screen('Close', [TexCircle TexBG]);
WaitSecs(duration);


end