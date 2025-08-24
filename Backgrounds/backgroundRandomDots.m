% Creates a random dot background

clc,clear,close all

imgSize = [3508, 3508]; % A3
outputDir = 'BOS_background';
mkdir(outputDir);

%%

% Parameters:
numDots = 100000;   % Number of dots
r = 2;              % Radius of dots


img_white = ones(imgSize);  
[xx, yy] = meshgrid(1:imgSize(1), 1:imgSize(2));
for i = 1:numDots
    cx = randi(imgSize(1));
    cy = randi(imgSize(2));

    % Circle mask
    xMin = max(1, cx - r); xMax = min(imgSize(1), cx + r);
    yMin = max(1, cy - r); yMax = min(imgSize(2), cy + r);
    xLocal = xx(yMin:yMax, xMin:xMax);
    yLocal = yy(yMin:yMax, xMin:xMax);
    mask = (xLocal - cx).^2 + (yLocal - cy).^2 <= r^2;

    regionB = img_white(yMin:yMax, xMin:xMax);
    regionB(mask) = min(regionB(mask), 1);
    img_white(yMin:yMax, xMin:xMax) = regionB;

end

imwrite(img_white, fullfile(outputDir, 'RandomDots_white_bg_dark_dots.png'));
