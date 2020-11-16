%%% In tutorial 2 we quantified and statistically tested recurrence
%%% in one dataset. In this tutorial we will perform a two level
%%% statistical analysis with four dummy participants to test for
%%% TGM recurrence across all participants.

% Load previously created brain time TGM structures
load tutorial2_output

%% Create four dummy participants with identical data
TGM_subj1.bt_TGM  = bt_TGMquant; % brain time quantified TGMs
TGM_subj1.bt_data = bt_data;     % brain time electrophysiological data
TGM_subj1.ct_TGM  = ct_TGMquant; % clock time quantified TGMs
TGM_subj1.ct_data = ct_data;     % clock time electrophysiological data

TGM_subj2.bt_TGM  = bt_TGMquant;
TGM_subj2.bt_data = bt_data;
TGM_subj2.ct_TGM  = ct_TGMquant;
TGM_subj2.ct_data = ct_data;

TGM_subj3.bt_TGM  = bt_TGMquant;
TGM_subj3.bt_data = bt_data;
TGM_subj3.ct_TGM  = ct_TGMquant;
TGM_subj3.ct_data = ct_data;

TGM_subj4.bt_TGM  = bt_TGMquant;
TGM_subj4.bt_data = bt_data;
TGM_subj4.ct_TGM  = ct_TGMquant;
TGM_subj4.ct_data = ct_data;

%% Brain Time warped data
% First level statistics (single subject level)
cfg.mvpacfg         = cfg_mv;          %input previous mvpa light config structure
cfg.figure          = 'no'; 
cfg.numperms1       = 5;               %number of permutations on the first level (per participant)
cfg.statsrange      = [1 20];          %range of tested recurrence rates
cfg.clabel          = clabel;

% Loop through all participants
for subj = 1:4 
data = eval(strcat('TGM_subj',num2str(subj),'.bt_data'));
TGM = eval(strcat('TGM_subj',num2str(subj),'.bt_TGM'));
[bt_stats1{subj}] = bt_TGMstatslevel1(cfg,data,TGM); %bt_TGMstatslevel2 requires one cell for each participant
end

% Apply second level statistics
cfg.numperms2      = 1000;                    %number of second level Monte Carlo permutations
cfg.multiplecorr   = 'fdr';                   %multiple correction option
[bt_pval] = bt_TGMstatslevel2(cfg,bt_stats1); %output matrix contains p-values and associated frequencies
title('Brain Time warped data results (SIGNIFICANT at simulated 8 Hz)')

%% Untransformed Clock Time data
% First level statistics (single subject level)
% Loop through all participants
for subj = 1:4 
data = eval(strcat('TGM_subj',num2str(subj),'.ct_data'));
TGM = eval(strcat('TGM_subj',num2str(subj),'.ct_TGM'));
[ct_stats1{subj}] = bt_TGMstatslevel1(cfg,data,TGM); %bt_TGMstatslevel2 requires one cell for each participant
end

% Apply second level statistics
[ct_pval] = bt_TGMstatslevel2(cfg,ct_stats1); %output matrix contains p-values and associated frequencies
title('Untransformed Clock Time data results (NOT SIGNIFICANT at simulated 8 Hz)')