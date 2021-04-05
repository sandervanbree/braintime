function setup_braintime
% Add all folders to path, and check for Signal Processing Toolbox,
% FieldTrip and MVPA Light.

% Add all folders and subfolders to path
braintime_path = fileparts(fileparts(mfilename('fullpath')));
try
    addpath(genpath(fullfile(braintime_path,'setup')));
    addpath(genpath(fullfile(braintime_path,'utilities')));
    addpath(genpath(fullfile(braintime_path,'clocktobrain')));
    addpath(genpath(fullfile(braintime_path,'tutorial')));
    addpath(genpath(fullfile(braintime_path,'periodicity')));
    addpath(genpath(fullfile(braintime_path,'datacheck')));
    addpath(genpath(fullfile(braintime_path,'warpingsource')));
    addpath(genpath(fullfile(braintime_path,'external')));
    addpath(genpath(fullfile(braintime_path,'topography')));
    addpath(genpath(fullfile(braintime_path,'dipolesimulation')));
catch
    error('Unable to add folders and subfolders to path. Please add all folders and subfolders to get started.');
end
    disp('All folders were successfully added to path.')

% Check for Signal Processing Toolbox
if exist('dtw') ~= 2
error('Unable to locate MATLAB Signal Processing Toolbox. Please install it by going to Add-Ons -> Get Add-Ons -> use Add-On Explorer to find the Signal Processing Toolbox');   
else
disp('MATLAB Signal Processing Toolbox is running...');    
end

% Check for FieldTrip
if exist('ft_freqanalysis') ~= 2
error('Unable to locate FieldTrip. Please add it to path (using ft_defaults) or download it here: https://www.fieldtriptoolbox.org/download/');
else
disp ('FieldTrip is up and running...');
end

% Check for MVPA Light
if exist('mv_classify_timextime') ~= 2
error('Unable to locate MVPA Light. Please set it up (using startup.m) or download it here: https://github.com/treder/MVPA-Light');   
else
disp('MVPA Light is up and running...');    
end

% If we reach here, all checks have passed
disp('Brain Time Toolbox setup successful - happy warping! The tutorials are a good place to start.');
end