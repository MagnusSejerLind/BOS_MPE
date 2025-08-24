%% === BOS Pattern Generator ===
clc; clear;

% === Settings ===
imgSize = [3508, 3508];
outputDir = 'BOS_patterns';
% mkdir(outputDir);

%% === (a) Random Dot Pattern ===
numDots = 50000*2;
% minRadius = 3;
% maxRadius = 10;

    r = 4;
    r=r/2


% Create 3 base images: black bg, white bg, white bg w/ color dots
% img_black = zeros(imgSize);         % Black background, white/gray dots
img_white = ones(imgSize);          % White background, black/gray dots
% img_color = ones([imgSize, 3]);     % White background, RGB color dots

[xx, yy] = meshgrid(1:imgSize(2), 1:imgSize(1));

for i = 1:numDots
    cx = randi(imgSize(2));
    cy = randi(imgSize(1));
    % r = randi([minRadius, maxRadius]);


    % Random gray intensity
    intensity = 1;

    % Circle mask
    xMin = max(1, cx - r); xMax = min(imgSize(2), cx + r);
    yMin = max(1, cy - r); yMax = min(imgSize(1), cy + r);
    xLocal = xx(yMin:yMax, xMin:xMax);
    yLocal = yy(yMin:yMax, xMin:xMax);
    mask = (xLocal - cx).^2 + (yLocal - cy).^2 <= r^2;

    % % Black background with bright dots
    % regionA = img_black(yMin:yMax, xMin:xMax);
    % regionA(mask) = max(regionA(mask), intensity);
    % img_black(yMin:yMax, xMin:xMax) = regionA;

    % White background with dark dots
    regionB = img_white(yMin:yMax, xMin:xMax);
    regionB(mask) = min(regionB(mask), 1 - intensity);
    img_white(yMin:yMax, xMin:xMax) = regionB;

    % % Colorful dots
    % color = rand(1, 3); % Random RGB color
    % for c = 1:3
    %     regionC = img_color(yMin:yMax, xMin:xMax, c);
    %     regionC(mask) = color(c);
    %     img_color(yMin:yMax, xMin:xMax, c) = regionC;
    % end
end

% imwrite(img_black, fullfile(outputDir, 'a1_black_bg_gray_dots.png'));
imwrite(img_white, fullfile(outputDir, 'a2_white_bg_dark_dots.png'));
% imwrite(img_color, fullfile(outputDir, 'a3_white_bg_colored_dots.png'));

%% === (b) Random Dot Pattern + Gaussian Noise Distortion ===
[X, Y] = meshgrid(1:imgSize(2), 1:imgSize(1));
scales = [2, 8, 32]; 
weights = [0.5, 0.3, 0.2];
dispFieldX = zeros(imgSize); dispFieldY = zeros(imgSize);

for i = 1:length(scales)
    sz = round(imgSize / scales(i));
    noiseX = imresize(randn(sz), imgSize, 'bilinear');
    noiseY = imresize(randn(sz), imgSize, 'bilinear');
    dispFieldX = dispFieldX + weights(i) * noiseX;
    dispFieldY = dispFieldY + weights(i) * noiseY;
end

% Stronger distortion
displacementFactor = 500;
u = displacementFactor * (mat2gray(dispFieldX) - 0.5);
v = displacementFactor * (mat2gray(dispFieldY) - 0.5);

% Apply warping
[Xq, Yq] = meshgrid(1:imgSize(2), 1:imgSize(1));
img_black_warped = interp2(X, Y, double(img_black), Xq + u, Yq + v, 'linear', 0);
imwrite(img_black_warped, fullfile(outputDir, 'b1_black_gaussian_warp.png'));

img_white_warped = interp2(X, Y, double(img_white), Xq + u, Yq + v, 'linear', 1);
imwrite(img_white_warped, fullfile(outputDir, 'b2_white_gaussian_warp.png'));

% For RGB image (colorful)
img_color_warped = zeros(size(img_color));
for c = 1:3
    img_color_warped(:,:,c) = interp2(X, Y, double(img_color(:,:,c)), Xq + u, Yq + v, 'linear', 1);
end
imwrite(img_color_warped, fullfile(outputDir, 'b3_color_gaussian_warp.png'));

%% === (c) Sinusoidal Fringe Pattern ===
[x, y] = meshgrid(1:imgSize(2), 1:imgSize(1));
freqU = 100; freqV = 100;
imgC = sin(2*pi*freqU*x/imgSize(2)) + sin(2*pi*freqV*y/imgSize(1));
imgC = mat2gray(imgC);
imwrite(imgC, fullfile(outputDir, 'pattern_c_sine_superposition.png'));

% %% === (d) Random Dot Distorted by Wavelet-Style Noise (Coiflet-based) ===
% % Choose wavelet type and level
% wname = 'db1'; % higher orders = smoother noise
% level = 100; % Reduce level to get finer structures
% waveletScale = 50; % Adjust to control warp strength (increase waveletScale to exaggerate distortion)
% 
% % Generate base random noise image
% noiseImg = rand(imgSize);
% 
% % Perform discrete wavelet decomposition
% [C, S] = wavedec2(noiseImg, level, wname);
% 
% % Zero out detail coefficients to keep low-frequency content only
% C_denoised = C;
% % Approximation coeffs at final level are untouched
% approx_len = prod(S(1,:));
% % Zero all detail coefficients
% C_denoised(approx_len+1:end) = 0;
% 
% % Reconstruct smooth wavelet noise
% waveletSmooth = waverec2(C_denoised, S, wname);
% waveletSmooth = mat2gray(waveletSmooth);
% 
% % Use its gradient as displacement field
% [gradX, gradY] = gradient(waveletSmooth);
% 
% u_wave = waveletScale * gradX;
% v_wave = waveletScale * gradY;
% 
% % Warp original random dot image
% imgD_B_wavelet = interp2(X, Y, double(img_black), Xq + u_wave, Yq + v_wave, 'linear', 0);
% % Warp original random dot image
% imgD_W_wavelet = interp2(X, Y, double(img_white), Xq + u_wave, Yq + v_wave, 'linear', 1);
% % Warp original random dot image
% imgD_C_wavelet = zeros(size(img_color));
% for c = 1:3
%     imgD_C_wavelet(:,:,c) = interp2(X, Y, double(img_color(:,:,c)), Xq + u_wave, Yq + v_wave, 'linear', 1);
% end
% % Display & save
% imshow(imgD_B_wavelet, []); title('(d) Random Dot + Wavelet Distortion (Black)');
% imwrite(imgD_B_wavelet, fullfile(outputDir, 'pattern_d_random_dot_wavelet_black.png'));
% % Display & save
% imshow(imgD_W_wavelet, []); title('(d) Random Dot + Wavelet Distortion (White)');
% imwrite(imgD_W_wavelet, fullfile(outputDir, 'pattern_d_random_dot_wavelet_white.png'));
% % Display & save
% imshow(imgD_C_wavelet, []); title('(d) Random Dot + Wavelet Distortion (Colour)');
% imwrite(imgD_C_wavelet, fullfile(outputDir, 'pattern_d_random_dot_wavelet_colour.png'));

%% === (f) Binary Fringe Pattern (BWFTP) ===

blockSize = 30;  % Adjust frequency via block size
rows = imgSize(1);
cols = imgSize(2);

% Create checkerboard by tiling
rowPattern = repmat([1 0], 1, ceil(cols/(2*blockSize)));
colPattern = repmat([1; 0], ceil(rows/(2*blockSize)), 1);
checkerPattern = kron(colPattern * rowPattern, ones(blockSize));

% Crop to image size
checkerPattern = checkerPattern(1:rows, 1:cols);

% Create inverse pattern (white squares on black)
inversePattern = 1 - checkerPattern;

% Save both
imwrite(checkerPattern, fullfile(outputDir, 'f1_checker_black_on_white.png'));
imwrite(inversePattern, fullfile(outputDir, 'f2_checker_white_on_black.png'));
