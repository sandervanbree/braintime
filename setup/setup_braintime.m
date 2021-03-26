function setup_braintime

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

% Add all folders and subfolders to path
braintime_path = fileparts(fileparts(mfilename('fullpath')));
try
    addpath(genpath(braintime_path));       
    addpath(genpath(fullfile(braintime_path,'setup')));
    addpath(genpath(fullfile(braintime_path,'clocktobrain')));
    addpath(genpath(fullfile(braintime_path,'tutorial')));
    addpath(genpath(fullfile(braintime_path,'TGMrecurrence')));
    addpath(genpath(fullfile(braintime_path,'datacheck')));
    addpath(genpath(fullfile(braintime_path,'warpingsource')));
    addpath(genpath(fullfile(braintime_path,'dependencies')));
    addpath(genpath(fullfile(braintime_path,'topography')));
catch
    error('Unable to add folders and subfolders to path. Please add all folders and subfolders to get started.');
end
    disp('All folders were successfully added to path.')


% If we reach here, all checks have passed
disp('Brain Time Toolbox setup successful - time to warp! The tutorials are a good place to start.');
end