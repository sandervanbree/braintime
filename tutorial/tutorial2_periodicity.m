%%% In tutorial 2 we will take the transformed brain time data and
%%% analyze it by applying a classifier, and cross-time generalization.
%%% We will analyze the periodicity in classifier performance to visually
%%% demonstrate that brain time warping has accentuated the simulated
%%% periodic pattern.

load tutorial1_output

% When testing for periodicity with the toolbox, it is important that the
% data is not high pass filtered (or done with extreme care), as this can
% introduce strong artefacts in patterns of classification performance.
% See van Driel et al., 2021; https://doi.org/10.1016/j.jneumeth.2021.109080
% We find that these patterns can be highly periodic, flawing the analysis.

% Timelock the data
cfg = [];
cfg.keeptrials     = 'yes';
cfg.removemean     = 'no';
ct_data            = ft_timelockanalysis(cfg,ct_data);
bt_data            = ft_timelockanalysis(cfg,bt_warpeddata.data);

% Extract classification labels (clabels)
clabel             = bt_warpeddata.clabel;

%% Create Clock and Brain Time TGMs
% Use MVPA Light to generate time generalization matrices (TGM)
cfg_mv.classifier  = 'lda';     % Linear Discriminant Analysis
cfg_mv.metric      = 'acc';     % Accuracy
cfg_mv.repeat      = 3;         % Number of repetitions; use higher number for real data
cfg_mv.cv          = 'kfold'; 
cfg_mv.k           = 5;         % Number of folds
[~, ct_mv]         = mv_classify_timextime(cfg_mv, ct_data.trial, clabel);
[~, bt_mv]         = mv_classify_timextime(cfg_mv, bt_data.trial, clabel); 

% Alternatively, you may be interested in classification without temporal
% generalization - i.e. what we abbreviate as the "diagonal" of TGMs. This
% is appropriate when the neural code is expected to evolve throughout
% trials. To analyze the diagonal using the toolbox, simply change the function
% called in MVPA Light (see tutorial 4 for more details):
%  [~, ct_mv]       = mv_classify_across_time(cfg_mv, ct_data.trial, clabel);
%  [~, bt_mv]       = mv_classify_across_time(cfg_mv, bt_data.trial, clabel); 

%% Quantify TGM periodicity (compare clock and brain time)
cfg = [];
cfg.bt_warpeddata   = bt_warpeddata;              % The warped data structure is required input
cfg.figure          = 'yes';
cfg.method          = 'tgm';                      % perform analysis over TGM ('tgm'), or 
                                                  % autocorrelation map of the TGM ('ac').
                                                  % Analyzing on the TGM is better when patterns
                                                  % are predicted in only parts of the TGM.
                                                  % Analyzing the AC map accentuates patterns
                                                  % present in most parts of the TGM.
                                                  
 cfg.periodicityfoi   = [1 22];                   % Range of tested periodicity rates

cfg.refdimension    = 'clocktime';                % Quantify periodicity as a function of seconds in the data
ct_quant         = bt_quantify(cfg,ct_mv);        % Clock time
title('Clock time periodicity');

cfg.refdimension    = 'braintime';                % Quantify periodicity as a function of passed cycles in the data
bt_quant         = bt_quantify(cfg,bt_mv);        % Brain time
title('Brain time periodicity');

% If you chose the right independent component, you should be able to see
% improvements in the autocorrelation map and periodicity spectrum after
% brain time warping.

% Save results for tutorial 3
save tutorial2_output
