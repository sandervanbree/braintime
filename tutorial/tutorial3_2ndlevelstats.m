%%% In tutorial 2 we quantified and statistically tested recurrence
%%% in one dataset. In this tutorial we will perform a two level
%%% statistical analysis with four dummy participants to test for
%%% TGM recurrence across all participants.

% Load previously created brain time TGM structures
load bt_TGM

%% Create four dummy participants with identical data
TGM_subj1.TGM     = bt_TGMquant;
TGM_subj1.bt_data = bt_data;

TGM_subj2.TGM     = bt_TGMquant;
TGM_subj2.bt_data = bt_data;

TGM_subj3.TGM     = bt_TGMquant;
TGM_subj3.bt_data = bt_data;

TGM_subj4.TGM     = bt_TGMquant;
TGM_subj4.bt_data = bt_data;

%% First level statistics (single subject level)
cfg.mvpacfg         = cfg_mv;          %input previous mvpa light config structure
cfg.figure          = 'no'; 
cfg.numperms1       = 5;               %number of permutations on the first level (per participant)
cfg.statsrange      = [1 20];          %range of tested recurrence rates
cfg.clabel          = clabel;

% Loop through all participants
for subj = 1:4 
data = eval(strcat('TGM_subj',num2str(subj),'.bt_data'));
TGM = eval(strcat('TGM_subj',num2str(subj),'.TGM'));
[stats1{subj}] = bt_TGMstatslevel1(cfg,data,TGM); %For bt_TGMstatslevel2 level to work, create one cell for each participants
end

%% Apply second level statistics
cfg.numperms2      = 100;
[pval] = bt_TGMstatslevel2(cfg,stats1); %Output matrix contains p-values and associated frequencies