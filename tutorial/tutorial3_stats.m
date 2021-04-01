%%% In tutorial 3 we will statistically analyze periodicity in the TGMs
%%% by performing a two-level statistical analysis with four dummy
%%% participants. Recurrence is compared against periodicity in
%%% data with permuted classification labels (serving as the null
%%% distribution).

% Load previously created brain time TGM structures
load tutorial2_output

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
cfg.mvpacfg         = cfg_mv;          % Input previous MVPA Light config structure
cfg.figure          = 'yes';           % Plot first level stats of each participant
cfg.numperms1       = 5;               % Number of permutations on the first level (per participant).
                                       % For real analyses, this should be higher.
cfg.statsrange      = [1 22];          % Range of tested periodicity rates
cfg.clabel          = clabel;          % We've saved clabel from last tutorial

for subj = 1:4 % This takes a minute or two for this data
data = eval(strcat('subj',num2str(subj),'.ct_data'));
quant = eval(strcat('subj',num2str(subj),'.ct_quant'));
[ct_stats1{subj}] = bt_statslevel1(cfg,data,quant); %bt_statslevel2 requires one cell for each participant
end

% The first level results will be very similar for all participants, as we
% just copy pasted datasets. Moreover, since our cfg.numperm1 is so low,
% significant testing is disabled on the first level. Let's just close 
% the plots and turn off 1st level plotting for now.
close all;
cfg.figure = 'no';

% Apply second level statistics (group-level)
cfg.numperms2      = 100000;                    % Number of second level Monte Carlo permutations
cfg.multiplecorr   = 'fdr';                     % Multiple correction option
cfg.cluster_p      = 0.05;                      % Threshold for TGM cluster significance testing
cfg.cluster_n      = 10;                        % Maximum number of clusters
cfg.cluster_smooth = 2;                         % Width of smoothing window (Gaussian SD), used only for cluster testing
[ct_stats2] = bt_statslevel2(cfg,ct_stats1);    % Output matrix contains p-values and associated frequencies
disp('Clock Time results (LOW periodicity)')

%% Brain time statistics
% If the right warping source was selected, brain time unveils simulated 
% periodicity. In addition, the time axis transforms from seconds to cycles,
% and the frequency space normalizes to participants' warping frequency
% (1 Hz = each participant's warping frequency). In addition, the toolbox
% displays significant clusters in the TGM.

for subj = 1:4 
data = eval(strcat('subj',num2str(subj),'.bt_data'));
quant = eval(strcat('subj',num2str(subj),'.bt_quant'));
[bt_stats1{subj}] = bt_statslevel1(cfg,data,quant);
end

[bt_stats2] = bt_statslevel2(cfg,bt_stats1); 
disp('Brain Time results (HIGH periodicity)')

%% Results extraction
disp(bt_stats2.warpedfreqs);   % Displays the warping frequency and its harmonics.
                               % Results are predicted at one of these frequencies.
                               
disp(bt_stats2.mapmethod);     % Was the analysis performed over the TGM, its
                               % AC map, or the diagonal?

disp(bt_stats2.periodicity);   % A vector of the second level periodicity spectra;
                               % their frequency and corrected p-values.
                               % Specifically the bt_stats2.warpedfreqs are
                               % exempted from multiple testing correction,
                               % because periodicity is predicted at those
                               % rates.
                               
disp(bt_stats2.clusters);      % These are the clusters in the tested dimension,
                               % and their associated p-values.

