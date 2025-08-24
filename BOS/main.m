% Background oriented schlieren method analysis
% Calculates the displacement and reflective index from provided data

clc,clear,close all
set(0,'defaultTextInterpreter','latex');

%% Pars

% Define the folder path
folderPath = 'Data/temp';

% BCs
opt.BC = 0; % [1/0]: Boundary condtions, 1: top/bottom, 0: left/right
   
% Pixel scaling value
scaleVal = 0.9;   

% interrogation window HALF size
wxy = 32;

% window spacing
sxy = 16;

%% Import data

% Import data
[lv_pathstr] = fileparts(which('readimx'));     % Define path to read function
lv_dirlist = dir([lv_pathstr '/' folderPath]);  % Find images to read

pix = struct();
lv_j = 1;
i = 0;
data = [];
for lv_i = 1:size(lv_dirlist,1)
    if (~lv_dirlist(lv_i).isdir) && (lv_dirlist(lv_i).bytes > 0) % check if item is valid
        i = i + 1;
        disp(['Loading ' lv_dirlist(lv_i).name]);
        currFile = readimx([lv_pathstr '/' folderPath '/' lv_dirlist(lv_i).name]);
        % figure; showimx(data(lv_j).Frames{1}); title(['Image ' num2str(lv_j) ': ' lv_dirlist(lv_i).name]);

        data= [data, currFile];
        pix(i).vals = currFile.Frames{1}.Components{1}.Planes{1};  % pixel vals
    end
end

% display first image (valid)
% figure()
% showimx(data(1).Frames{1}); title(['First Valid Image: ' lv_dirlist(find([lv_dirlist.bytes] > 0, 1)).name]);

%% encoding and scaling

for i = 1:length(pix)
    img = pix(i).vals;

    % Normalize and convert to uint16
    img = img - min(img(:));
    img = img / (max(img(:)) / scaleVal );
    img = uint16(img * 65535 );
    img = imrotate(img, -90);

    pix(i).vals = img;
end

%% define analysis region

figure()
disp('Select region of interest with cursor - And double click region');
[as,rect]= imcrop(pix(1).vals);
bs= imcrop(pix(end).vals,rect);


%% Calculate displacements
[x,y,uwarp,vwarp,iter]= FindDisplacements (as,bs,[wxy,wxy],[sxy,sxy]);

mag = sqrt(uwarp.^2 + vwarp.^2);    % Magnitude of displacements
% mag(mag > 5) = nan;   % apply to remove glare
%% Plot displacement

% quiver plot
figure()
quiver(x,y,vwarp,uwarp,'k','AutoscaleFactor', 3);
title('Displacement')
xlabel('x [mm]')
ylabel('y [mm]')
xlim([min(x(:)) max(x(:))])
ylim([min(y(:)) max(y(:))])

% quiver plot - fewer data points included
figure()
step = 2; 
quiver(x(1:step:end,1:step:end), y(1:step:end,1:step:end), vwarp(1:step:end,1:step:end), uwarp(1:step:end,1:step:end), 'k', 'AutoscaleFactor', 5);
title('Displacement')
xlabel('x [mm]')
ylabel('y [mm]')
xlim([min(x(:)) max(x(:))])
ylim([min(y(:)) max(y(:))])


% scatter plot
figure()
axis ij
scatter(x(:), y(:), 15, mag(:), 'filled');
colormap(jet)
colorbar
title('Displacement')
xlabel('x [mm]')
ylabel('y [mm]')
xlim([min(x(:)) max(x(:))])
ylim([min(y(:)) max(y(:))])
% clim([0, 1]);


% contour plot
figure()
contourf(x, y, mag)
colormap(jet)
colorbar
title('Displacement')
% clim([0, 1]);


%% Determine refractive index

% fluid properties - air
bri = 1.000273;     % Base Refractive Index STP
bri_dT = -0.0001;   % Rate of Change of Refractive Index 


% BC - air at STP
BC.A = bri*ones(size(uwarp(1,:)));
BC.B = bri*ones(size(uwarp(end,:)));

% BC - air non-STP
    % temp = 40; % [°c]
    % bc_n = (temp-20)*bri_dT + bri;
    % BC.A = bc_n*ones(size(uwarp(:,end)));
    % BC.B = BC.A;


% integrate the refractive index gradient field
nxy = ScaleAndIntegrate(uwarp,vwarp,BC,opt);

%% Plot reflective index

figure()
contourf(x, y, nxy, 25,'LineColor', 'none');
colormap jet
colorbar
title ('Refractive Index Distribution')
xlabel('x [mm]')
ylabel('y [mm]')
xlim([min(x(:)) max(x(:))])
ylim([min(y(:)) max(y(:))])
clim([0.9993, 1.0005])
