function setup_braintime
%Automatically add all subfolders to path
braintime_path = fileparts(fileparts(mfilename('fullpath')));

addpath(braintime_path);
addpath(fullfile(braintime_path,'setup'));
addpath(fullfile(braintime_path,'clocktobrain'));
addpath(fullfile(braintime_path,'tutorial'));
addpath(fullfile(braintime_path,'TGMrecurrence'));
addpath(fullfile(braintime_path,'datacheck'));
addpath(fullfile(braintime_path,'warpingsource'));
addpath(fullfile(braintime_path,'dependencies'));
addpath(fullfile(braintime_path,'topography'));

disp('All folders added to path successfully. See tutorials to get started.')