%%% In tutorial 4 we will show that besides testing periodicity in TGMs,
%%% you can also test for periodicity on the diagonal of TGMs. That is,
%%% analyzing classification when training and testing on the same
%%% timepoint, ignoring how training on one timepoint generalizes to other
%%% timepoints. Tutorial 4 differs minimally from tutorial 2 and 3 -- the
%%% only difference is the analysis is done over a 1D vector (the diagonal
%%% of the TGM), rather than a 2D matrix (the TGM). Theoretically, this
%%% analysis is optimal when the neural signature underlying your classes
%%% evolves within your time window of interest.

% Let's start off with the variables we had right before bt_quantify
load tutorial2_output cfg_mv ct_data bt_data clabel bt_warpeddata

% mv_classify_timextime (TGM) ---> mv_classify_across_time (Diagonal)
[~, ct_mv]       = mv_classify_across_time(cfg_mv, ct_data.trial, clabel);
[~, bt_mv]       = mv_classify_across_time(cfg_mv, bt_data.trial, clabel); 

%% Quantify TGM periodicity (compare clock and brain time)
cfg = [];
cfg.bt_warpeddata   = bt_warpeddata;              % The warped data structure is required input
cfg.figure          = 'yes';
% cfg.method          = 'tgm';                    % This is not needed - there is only one method for
                                                  % diagonal analysis
cfg.periodicityfoi   = [1 22];                    % Range of tested periodicity rates

cfg.refdimension    = 'clocktime';                % Quantify periodicity as a function of seconds in the data
ct_quant            = bt_quantify(cfg,ct_mv);    % Clock time
title('Clock time periodicity');

cfg.refdimension    = 'braintime';                % Quantify periodicity as a function of passed cycles in the data
bt_quant            = bt_quantify(cfg,bt_mv);    % Brain time
title('Brain time periodicity');

% Note how recurrence may be higher than the warped frequency. Depending
% on the underlying process, periodicity may be expected at 0.5x, 1x or 2x
% the warped frequency. See the Brain Time Toolbox paper for more details.

%% Create four dummy participants with identical data
subj1.bt_quant  = bt_quant; % brain time quantified TGMs
subj1.bt_data   = bt_data;  % brain time electrophysiological data
subj1.ct_quant  = ct_quant; % clock time quantified TGMs
subj1.ct_data   = ct_data;  % clock time electrophysiological data

subj2.bt_quant  = bt_quant;
subj2.bt_data   = bt_data;
subj2.ct_quant  = ct_quant;
subj2.ct_data   = ct_data;

subj3.bt_quant  = bt_quant;
subj3.bt_data   = bt_data;
subj3.ct_quant  = ct_quant;
subj3.ct_data   = ct_data;

subj4.bt_quant  = bt_quant;
subj4.bt_data   = bt_data;
subj4.ct_quant  = ct_quant;
subj4.ct_data   = ct_data;

%% Clock time statistics
% Let's first show that the original data shows little to no periodicity

% First level statistics (single subjects)
cfg.figure          = 'no'; 
cfg.numperms1       = 5;               % Number of permutations on the first level (per participant).
                                       % For real analyses, this should be higher.
cfg.statsrange      = [1 22];          % Range of tested periodicity rates

for subj = 1:4 % This takes a minute or two for this data
data = eval(strcat('subj',num2str(subj),'.ct_data'));
quant = eval(strcat('subj',num2str(subj),'.ct_quant'));
[ct_stats1{subj}] = bt_statslevel1(cfg,data,quant); 
end

% Apply second level statistics
cfg.numperms2      = 100000;                    % Number of second level Monte Carlo permutations
cfg.multiplecorr   = 'fdr';                     % Multiple correction option

% By the way, many cfg_clus parameters are set to default values when not set manually
cfg_clus.test             = 'permutation';            
cfg_clus.correctm         = 'cluster';           
cfg_clus.n_permutations   = 1000;               
[ct_stats2] = bt_statslevel2(cfg,cfg_clus,ct_stats1);    
disp('Clock Time results (LOW periodicity)')

%% Brain time statistics
for subj = 1:4 
data = eval(strcat('subj',num2str(subj),'.bt_data'));
quant = eval(strcat('subj',num2str(subj),'.bt_quant'));
[bt_stats1{subj}] = bt_statslevel1(cfg,data,quant);
end

[bt_stats2] = bt_statslevel2(cfg,cfg_clus,bt_stats1); 
disp('Brain Time results (HIGH periodicity)')

% Depending on the data structure, the TGM and diagonal may show recurrence
% not at 1 Hz, but at half (0.5 Hz) or double (2 Hz) the warping frequency.

% In this simulation, spurious significant frequencies may arise as one
% small dataset is copied four times over.