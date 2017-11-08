%This function shows fixation circles of given color at a selected position on the screen.
%
% Input:
% win     = window Pointer
% winrect = window rectangle, as returned from OpenWindow
% sizes   = vector of circle-diameters. The drawing follows the order in
% the vector, so it is required to start with the largest one
% colors  = vector of circle-colors (in greyscale) or matrix of RGB values
% with every column being one circle rgb
%
%   Example showing three circles of size 10,6,2
%
%       whichScreen = max(Screen('Screens'));
%       [win winRect] = Screen('OpenWindow',whichScreen,0);
%       Screen('FillRect', win, [128 128 128 1]);
%       ShowFixationCircles(win,winRect,[10,6,2],[0,128,0])
%       Screen('Flip', win);
%

function ShowFixationCircles(win,winRect,sizes,colors,x,y)

%if not specified, show it in the center
if nargin <5
    x = winRect(3)/2;
    y = winRect(4)/2;
end

%calculate the bounding boxes and set the colors
for k=1:length(sizes)
    rect(:,k)=[x-sizes(k)/2 y-sizes(k)/2 x+sizes(k)/2 y+sizes(k)/2];
end

if size(colors,1)==1;
    fx_colors = repmat(colors,3,1);
else
    fx_colors = colors;
end

Screen( 'FillOval', win, fx_colors, rect);

end
