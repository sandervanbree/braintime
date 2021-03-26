function setup_braintime
%Automatically add all subfolders to path
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
    disp('(Error while adding folders to path. Please add all folders and subfolders to get started.');
    disp('Then, see tutorials to get started.');
end
    disp('All folders added to path successfully. See tutorials to get started.')
end