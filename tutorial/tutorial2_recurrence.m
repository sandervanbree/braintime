%%% In tutorial 2 we will take the transformed brain time data and
%%% analyze recurrence present in it by applying cross-time
%%% generalization.

% Load previously created clock and brain time data
load tutorial1_output

% Timelock analyze the data
cfg = [];
cfg.keeptrials     = 'yes';
cfg.removemean     = 'no';
ct_data            = ft_timelockanalysis(cfg,ct_data);
bt_data            = ft_timelockanalysis(cfg,bt_warpeddata.data);

% Extract classification labels (clabels)
clabel             = bt_warpeddata.clabel;

%% Create Clock and Brain Time TGM
% Use MVPA Light to visualize recurrence in a time generalization matrix (TGM)
cfg_mv.classifier  = 'lda';
cfg_mv.metric      = 'acc';     % Accuracy
cfg_mv.repeat      = 3;         % Number of repetitions; use higher than 1 for real data
cfg_mv.cv          = 'kfold';
cfg_mv.k           = 5;         % Number of folds
[ct_TGM, ~]        = mv_classify_timextime(cfg_mv, ct_data.trial, clabel);
[bt_TGM, ~]        = mv_classify_timextime(cfg_mv, bt_data.trial, clabel); 

%% Quantify TGM recurrence (compare clock and brain time)
cfg = [];
cfg.bt_warpeddata   = bt_warpeddata;              % Specify so that information can be retrieved
cfg.MVPAcfg         = cfg_mv;                     % Input MVPA light config structure
cfg.figure          = 'yes';
cfg.mapmethod       = 'ac';                       % perform analysis over TGM's autocorrelation map
cfg.recurrencefoi   = [6 20];                     % Range of tested recurrence rates

cfg.refdimension    = 'clocktime';                % Quantify recurrence as a function of seconds in the data
ct_TGMquant         = bt_TGMquantify(cfg,ct_TGM); % Compare with clock time
title('Clock time recurrence');

cfg.refdimension    = 'braintime';                % Quantify recurrence as a function the warped frequency
bt_TGMquant         = bt_TGMquantify(cfg,bt_TGM); % Do once for brain time
title('Brain time recurrence');

%% Statistically test TGM recurrence on the single subject level (compare clock and brain time)
clabel = bt_warpeddata.clabel;

cfg.numperms1       = 5;               % Number of permutations on the first level
cfg.clabel          = clabel;
[ct_TGMstats1] = bt_TGMstatslevel1(cfg,ct_data,ct_TGMquant);  % Clock time results (not significant)
[bt_TGMstats1] = bt_TGMstatslevel1(cfg,bt_data,bt_TGMquant);  % Brain time results (significant)

%% Save results for tutorial 3
save tutorial2_output ct_TGMquant bt_TGMquant ct_TGMstats1 bt_TGMstats1 ct_data bt_data clabel cfg_mv
