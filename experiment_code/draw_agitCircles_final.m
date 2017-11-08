
function flipTimestamp = draw_agitCircles_final(screenRes, window, windowRect, targetVisAngle, px_per_deg, duration, number_circles, ctrDot, gaussCircles, gaussMask) %targetVisAngle [degrees], duration [s], ctrDot[boolean], gauss[boolean]

screenFR = Screen('NominalFrameRate', window);
[window_width, window_height] = Screen('WindowSize', window);

time1 = GetSecs();
time2 = GetSecs();
times = [];

%%Test velocity
% test = [];
% flipList = [];

velo = 3;   %velocity of stimulus refreshrate -> 1/velo * screenFR --> the higher the slower
s_per_flip = 1/screenFR;
timeBetween = s_per_flip * velo;
flipTimestamp = clock;

diamAngle = targetVisAngle; %diameter of target in degrees

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

phase = 0; %set phase to be changed in while loop for a "moving" sinus for obtaining the moving circles

%%calc sine values as gray indices
sin_upBound = pi * number_circles - pi/2;
tmp_idx = find(circleRect(:,:) == 1);

%Draw background
TexBG = Screen('MakeTexture', window, bigGrid);
Screen('DrawTexture', window, TexBG);

while((time2 - time1) < duration)
    tic;

    %get y-value of sinus-function of the respective distance to circle-center as x-value
    distance = (sqrt((X(tmp_idx)-ctr_XY).^2 + (Y(tmp_idx)-ctr_XY).^2)); %distance of each arrayfield from center
    circleRect(tmp_idx) = sin(distance*(sin_upBound/radiusCross) + phase);
    
    %%gaussian circles --> sinusvalues: 1=black & -1=white  
    if(gaussCircles == 1) 
      circleRect(tmp_idx) = (circleRect(tmp_idx) + 1).* colorFac; 
    else    %sharp circles
      circleRect(find(circleRect <= 0 & circleRect >= -1)) = white;
      circleRect(find(circleRect > 0 & circleRect <= 1)) = black; 
    end

% %     convParam = 0.2;
% %     matrix = [convParam convParam convParam; convParam 1 convParam; convParam convParam convParam];
%     matrix = matrix./sum(matrix(:));

%     matrix = fspecial('gaussian',11,1.5);
%     circleRect = background-filter2(matrix, background-circleRect,'same');
    
     %Draw Circle 
    TexCircle = Screen('MakeTexture', window, circleRect);
    Screen('DrawTexture', window, TexCircle, [], destRect);
        
    %%Gauss mask
    if(gaussMask == 1) %check if gauss argument is true
        mask = ones(rectSide, rectSide, 2);
        [Xm, Ym] = meshgrid(-ctr_XY+1:ctr_XY-1, -ctr_XY+1:ctr_XY-1);
        ctr_XYm = ctr_XY/6;
        mask(:,:,1) = background;
        mask(:,:,2) = round(255 - exp(-((Xm/ctr_XYm).^2)-((Ym/ctr_XYm).^2)).*255);
        mask2 = mask(:,:,2);
        TexMask = Screen('MakeTexture', window, mask);
        destRectMask = [0 0 rectSide rectSide];
        [destRectMask, dhm, dvm] = CenterRect(destRectMask, windowRect);
        Screen('DrawTexture', window, TexMask, [], destRectMask);
    else %smooth outer edge of circle
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
    
    %%Center Dot
        if(ctrDot == 1) %if ctrDot == true
            dotR = 2; % [px]
            dotHaloR = dotR*2;
             
            centerDot = ones(dotHaloR*2+1, dotHaloR*2+1)*background;
            [Xd Yd] = meshgrid(1: dotHaloR*2+1, 1:dotHaloR*2+1);
            ctr_centerDot = (dotHaloR*2+2)/2;
            
            centerDot(sqrt((Xd-ctr_centerDot).^2 + (Yd-ctr_centerDot).^2) <= dotR) = black;
            centerDot2 = centerDot*0;
            centerDot2(sqrt((Xd-ctr_centerDot).^2 + (Yd-ctr_centerDot).^2) <= dotHaloR) = 255; 
            centerDot = cat(3, centerDot, centerDot2);

            [cntrDotRect, dh, dv] = CenterRect([0 0 dotHaloR+1 dotHaloR+1], windowRect);
            
            centerDot = Screen('MakeTexture', window, centerDot);
            Screen('DrawTexture', window, centerDot, [], cntrDotRect);
        end

    times = [times toc];
    whenFlip = timeBetween - times(end); %get time between start of this while-turn and now and flip only at next timepoint when flip should happen (= every 1/timebetween)
    [vbl stimOn flipTimestamp] = Screen('Flip', window, (flipTimestamp(1) + whenFlip),1);
    %     test = [test (flipTimestamp - vbl)];
    %     flipList = [flipList flipTimestamp];
    Screen('Close', [TexMask TexBG  TexMask]);
    phase = phase + 1;
    time2 = GetSecs();
    
end


end





