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
bt_data            = ft_timelockanalysis(cfg,bt_struc.data);

% Extract classification labels (clabels)
clabel             = bt_struc.clabel;

%% Create Clock and Brain Time TGM
% Use MVPA Light to visualize recurrence in a time generalization matrix (TGM)
cfg_mv.classifier  = 'lda';
cfg.preprocess     = 'yes';     %z-scoring the data can improve classification
cfg_mv.metric      = 'acc';     %accuracy
cfg_mv.repeat      = 1;         %number of repetitions; use higher than 1 for real data
cfg_mv.cv          = 'kfold';
cfg_mv.k           = 5;         %number of folds
[ct_TGM, ~]        = mv_classify_timextime(cfg_mv, ct_data.trial, clabel);
[bt_TGM, ~]        = mv_classify_timextime(cfg_mv, bt_data.trial, clabel); 

% Optional: Smooth and Z-score TGMs
% ct_TGM = zscore(ct_TGM,0,'all');  
% ct_TGM = imgaussfilt(ct_TGM,0.2); %Smoothing TGM may accentuate recurrence
% bt_TGM = zscore(bt_TGM,0,'all');
% bt_TGM = imgaussfilt(bt_TGM,0.2);

% Plot results
figure; subplot(1,2,1)
mv_plot_2D(ct_TGM);title('Clock time TGM')
subplot(1,2,2)
mv_plot_2D(bt_TGM);title('Brain time TGM')

%% Quantify TGM recurrence (compare clock and brain time)
cfg = [];
cfg.bt_struc        = bt_struc;                   %specify so that information can be retrieved
cfg.figure          = 'yes';
cfg.refdimension    = 'braintime';                %quantify recurrence as a function of seconds in the data
ct_TGMquant         = bt_TGMquantify(cfg,ct_TGM); %compare with clock time
cfg.refdimension    = 'clocktime';                %quantify recurrence as a function the warped frequency
bt_TGMquant         = bt_TGMquantify(cfg,bt_TGM); %do once for brain time

%% Statistically test TGM recurrence on the single subject level (compare clock and brain time)
clabel = bt_struc.clabel;

cfg.mvpacfg         = cfg_mv;          %input previous mvpa light config structure
cfg.numperms1       = 5;              %number of permutations on the first level
cfg.statsrange      = [1 20];          %range of tested recurrence rates
cfg.clabel          = clabel;
[ct_TGMstats1] = bt_TGMstatslevel1(cfg,ct_data,ct_TGMquant);  %clock time results (not significant)
[bt_TGMstats1] = bt_TGMstatslevel1(cfg,bt_data,bt_TGMquant);  %brain time results (significant)


%% Save results for tutorial 3
save tutorial2_output ct_TGMquant bt_TGMquant ct_TGMstats1 bt_TGMstats1 ct_data bt_data clabel cfg_mv
