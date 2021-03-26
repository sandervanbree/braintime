%%% In tutorial 3 we will statistically analyze recurrence in the TGMs
%%% by performing a two-level statistical analysis with four dummy
%%% participants. Recurrence is compared against recurrence in
%%% data with permuted classification labels (serving as the null
%%% distribution).

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

%% Clock time statistics
% Let's first show that the original data shows little to no recurrence

% First level statistics (single subjects)
cfg.mvpacfg         = cfg_mv;          % Input previous MVPA Light config structure
cfg.figure          = 'yes';           % Plot first level stats of each participant
cfg.numperms1       = 5;               % Number of permutations on the first level (per participant).
                                       % For real analyses, this should be higher.
cfg.statsrange      = [1 22];          % Range of tested recurrence rates
cfg.clabel          = clabel;          % We've saved clabel from last tutorial

for subj = 1:4 % This takes a minute or two for this data
data = eval(strcat('TGM_subj',num2str(subj),'.ct_data'));
TGM = eval(strcat('TGM_subj',num2str(subj),'.ct_TGM'));
[ct_stats1{subj}] = bt_TGMstatslevel1(cfg,data,TGM); %bt_TGMstatslevel2 requires one cell for each participant
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
[ct_stats2] = bt_TGMstatslevel2(cfg,ct_stats1); % Output matrix contains p-values and associated frequencies
disp('Clock Time results (LOW recurrence)')

%% Brain time statistics
% If the right warping source was selected, brain time unveils simulated 
% recurrence. In addition, the time axis transforms from seconds to cycles,
% and the frequency space normalizes to participants' warping frequency
% (1 Hz = each participant's warping frequency). In addition, the toolbox
% displays significant clusters in the TGM.

for subj = 1:4 
data = eval(strcat('TGM_subj',num2str(subj),'.bt_data'));
TGM = eval(strcat('TGM_subj',num2str(subj),'.bt_TGM'));
[bt_stats1{subj}] = bt_TGMstatslevel1(cfg,data,TGM);
end

[bt_stats2] = bt_TGMstatslevel2(cfg,bt_stats1); 
disp('Brain Time results (HIGH recurrence)')

%% Results extraction
disp(bt_stats2.warpedfreqs);   % Displays the warping frequency and its harmonics.
                               % Results are predicted at one of these frequencies.
                               
disp(bt_stats2.mapmethod);     % Over what dimension was the analysis performed?

disp(bt_stats2.recurrence);    % A vector of the second level recurrence spectra;
                               % their frequency and corrected p-values.
                               % Specifically the bt_stats2.warpedfreqs are
                               % exempted from multiple testing correction,
                               % because recurrence is predicted at those
                               % rates.
                               
disp(bt_stats2.TGMclusters);   % These are the clusters in the tested dimension,
                               % and their associated p-values.

