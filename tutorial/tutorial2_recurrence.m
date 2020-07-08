%%% In tutorial 2 we will take the transformed brain time data and
%%% analyze recurrence present in it by applying cross-time
%%% generalization.

% Load previously created clock and brain time data (see tutorial folder)
load ct_data
load bt_struc

% Timelock analyze the data
cfg = [];
cfg.keeptrials     = 'yes';
cfg.removemean     = 'no';
ct_data            = ft_timelockanalysis(cfg,ct_data);
bt_data            = ft_timelockanalysis(cfg,bt_struc.data);

% Extract classification labels (clabels)
clabel             = bt_struc.clabel;

%% Create Clock and Brain Time TGM
% Use MVPA Light to visualize recurrence in a time generalization matrix (TGM)
cfg_mv.classifier  = 'lda';
cfg_mv.metric      = 'acc';     %accuracy
cfg_mv.repeat      = 3;         %number of repetitions
cfg_mv.cv          = 'kfold';
cfg_mv.k           = 5;         %number of folds
[ct_TGM, ~]        = mv_classify_timextime(cfg_mv, ct_data.trial, clabel);
[bt_TGM, ~]        = mv_classify_timextime(cfg_mv, bt_data.trial, clabel); 

% Optional: Smooth and Z-score TGMs
% ct_TGM = zscore(ct_TGM,0,'all');
% ct_TGM = imgaussfilt(ct_TGM,0.2);
% bt_TGM = zscore(bt_TGM,0,'all');
% bt_TGM = imgaussfilt(bt_TGM,0.2);

% Plot results
figure; subplot(1,2,1)
mv_plot_2D(bt_TGM);title('Brain time TGM')
subplot(1,2,2)
mv_plot_2D(ct_TGM);title('Clock time TGM')

%% Quantify TGM recurrence (compare clock and brain time)
cfg = [];
cfg.bt_struc        = bt_struc;      %specify so that information can be retrieved
cfg.refdimension    = 'braintime';   %quantify recurrence as a function of seconds in the data or the warped frequency
cfg.figure          = 'yes';
bt_quantTGM         = bt_quantifyTGM(cfg,bt_TGM); %do once for brain time
ct_quantTGM         = bt_quantifyTGM(cfg,ct_TGM); %compare with clock time

%% Statistically test TGM recurrence (compare clock and brain time)
cfg.mvpacfg         = cfg_mv;          %input previous mvpa light config structure
cfg.permlevels      = 1;               %for data with multiple participants, two level.
cfg.numperms1       = 10;              %number of permutations on the first level
cfg.statsrange      = [1 20];          %range of tested recurrence rates
cfg.clabel          = bt_struc.clabel;
bt_statsTGM(cfg,bt_data,bt_quantTGM);  %brain time results (significant)
bt_statsTGM(cfg,ct_data,ct_quantTGM);  %clock time results (not significant)

