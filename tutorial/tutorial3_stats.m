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
cfg.figure          = 'yes';           % Plot first level stats of each participant
cfg.numperms1       = 5;               % Number of permutations on the first level (per participant).
                                       % For real analyses, this should be higher.
cfg.statsrange      = [1 22];          % Range of tested periodicity rates

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
                                                % We recommend at least 100000, or ideally an order higher
cfg.multiplecorr   = 'fdr';                     % Multiple correction option ('none', 'fdr', 'bonferroni')

% A plethora of cluster parameters are available that are sent to MVPA
% Light for cluster correction. For details on each parameter, check out
% .../MVPA-Light-master/examples/understanding_statistics.m

cfg_clus.test             = 'permutation';  
cfg_clus.correctm         = 'cluster';           
cfg_clus.n_permutations   = 1000;                % Number of permutations
cfg_clus.clusterstatistic = 'maxsum';            % Sums cluster statistics
cfg_clus.tail             = 1;                   % One-sided test for above chance performance
cfg_clus.alpha            = 0.05;                % Alpha value
cfg_clus.statistic        = 'wilcoxon';          % Non-parametric counterpart to t-stats
cfg_clus.null             = 0.5;                 % By default, most performance
                                                 % metric's null is at 0.5
cfg_clus.clustercritval   = 1.96;                % The critical z-value 
% (From MVPA Light documentation):
% z-val = 1.65 corresponds to uncorrected p-value = 0.1
% z-val = 1.96 corresponds to uncorrected p-value = 0.05
% z-val = 2.58 corresponds to uncorrected p-value = 0.01

[ct_stats2] = bt_statslevel2(cfg,cfg_clus,ct_stats1);    % Output matrix contains p-values and associated frequencies
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

[bt_stats2] = bt_statslevel2(cfg,cfg_clus,bt_stats1); 
disp('Brain Time results (HIGH periodicity)')

%% Results extraction
disp(bt_stats2.warpedfreqs);   % Displays the warping frequency and its harmonics.
                               % Results are predicted at one of these frequencies.
                               
disp(bt_stats2.periodicity);   % A vector of the second level periodicity spectra;
                               % their frequency and corrected p-values.
                               % Specifically the bt_stats2.warpedfreqs are
                               % exempted from multiple testing correction,
                               % because periodicity is predicted at those
                               % rates.
                               
disp(bt_stats2.analysistype);  % Was the analysis performed over the TGM, its
                               % AC map, or the diagonal?

disp(bt_stats2.clusters);      % MVPA Light cluster results structure, with p-values
                               % and cluster indices
