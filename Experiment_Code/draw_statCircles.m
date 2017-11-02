
function draw_statCircles(screenRes, window, windowRect, targetVisAngle, px_per_deg, duration, number_circles, ctrDot, gaussMask) %targetVisAngle [degrees], duration [s], ctrDot[boolean], gauss[boolean]

screenFR = Screen('NominalFrameRate', window);
[window_width, window_height] = Screen('WindowSize', window);

% time1 = GetSecs();
% time2 = GetSecs();
% times = [];

diamAngle = targetVisAngle; %diameter of target in degrees: for a subject-monitor distance of 600mm!!

diamCross = diamAngle * px_per_deg;
radiusCross = round(0.5*diamCross);
radiusDelta = round(radiusCross/number_circles); %has to conform with the number of defined radiusses (see below)

%%Set black/white & gray
white = WhiteIndex(0);
black = BlackIndex(0);
colorFac = white/2;
background = 1*colorFac;

bigGrid(1:window_height,1:window_width) = background;
ctr_Grid = size(bigGrid)/2;
puffer = 80; %necessary for gauss mask

rectSide = 2*radiusCross+puffer+1; %size of side of the rectangle the gauss mask will be applied to
circleRect = ones(rectSide, rectSide)*background; %Rectangle to plot circles in

destRect = [0 0 (rectSide) (rectSide)];
[destRect, dh, dv] = CenterRect(destRect, windowRect);

%Define which fields are lying within the circle
X = repmat((1:rectSide), rectSide, 1);
Y = repmat((1:rectSide)', 1, rectSide);
ctr_XY = (rectSide+1)/2;
circleRect(sqrt((X-ctr_XY).^2 + (Y-ctr_XY).^2) <= radiusCross) = 1;

%%calc sine values as gray indices
sin_upBound = pi * number_circles - pi/2;
tmp_idx = find(circleRect(:,:) == 1);

%Draw background
TexBG = Screen('MakeTexture', window, bigGrid);
Screen('DrawTexture', window, TexBG);

%get y-value of sinus-function of the respective distance to circle-center as x-value
distance = (sqrt((X(tmp_idx)-ctr_XY).^2 + (Y(tmp_idx)-ctr_XY).^2));
circleRect(tmp_idx) = sin(distance*(sin_upBound/radiusCross));

%%gaussian circles --> sinusvalues: 1=black & -1=white
%       circleRect(tmp_idx) = (circleRect(tmp_idx) + 1).* colorFac;

%%sharp circles
circleRect(find(circleRect <= 0 & circleRect >= -1)) = white;
circleRect(find(circleRect > 0 & circleRect <= 1)) = black;

%Draw Circle
TexCircle = Screen('MakeTexture', window, circleRect);
Screen('DrawTexture', window, TexCircle, [], destRect);
        
    %%Gauss mask (for sharp circles)
    if(gaussMask == 1) %check if gauss argument is true
        mask = ones(rectSide, rectSide, 2);
        [Xm, Ym] = meshgrid(-ctr_XY+1:ctr_XY-1, -ctr_XY+1:ctr_XY-1);
        ctr_XYm = ctr_XY/6;
        mask(:,:,1) = background;
        mask(:,:,2) = round(255 - exp(-((Xm/ctr_XYm).^2)-((Ym/ctr_XYm).^2)).*255);
        TexMask = Screen('MakeTexture', window, mask);
        destRectMask = [0 0 rectSide rectSide];
        [destRectMask, dhm, dvm] = CenterRect(destRectMask, windowRect);
        Screen('DrawTexture', window, TexMask, [], destRectMask);
    else    %smooth outer edge of circle
        mask1 = ones(rectSide, rectSide)*background;
        mask2 = ones(rectSide, rectSide);
        mask2(sqrt((X-ctr_XY).^2 + (Y-ctr_XY).^2) <= radiusCross-5) = 0;
        mask2(sqrt((X-ctr_XY).^2 + (Y-ctr_XY).^2) >= radiusCross) = 255;
        tmp_idx2 = find(mask2(:,:) == 1);
        distance2 = (sqrt((X(tmp_idx2)-ctr_XY).^2 + (Y(tmp_idx2)-ctr_XY).^2));
        mask2(tmp_idx2) = (distance2 - min(distance2)) .* (255/(max(distance2)-min(distance2)));
        mask = cat(3, mask1, mask2);
        TexMask = Screen('MakeTexture', window, mask);
        destRectMask = [0 0 rectSide rectSide];
        [destRectMask, dhm, dvm] = CenterRect(destRectMask, windowRect);
        Screen('DrawTexture', window, TexMask, [], destRectMask);
    end
    
 
Screen('Flip', window);
Screen('Close', [TexCircle TexBG TexMask]);
WaitSecs(duration);


end





