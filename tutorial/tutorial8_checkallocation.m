%%% In tutorial 7, we will compare the two cutting
%%% methods ('cutartefact' and 'consistenttime'), and see whether they
%%% put the same data in the same cycles. This is an important data control
%%% procedure, for users who want to compare the effect of the first cycle
%%% artefact in the TGM. While cutartefact cuts this artefact out, the
%%% longer time window implemented may destroy the data-to-cycle allocation
%%% that the toolbox applied in consistenttime.

%% Step 1: Load and structure data
load dipolesim_tutorial;

cfg               = [];
ct_data           = ft_appenddata(cfg,ct_left,ct_right);
clabel            = [ones(size(ct_left.trial,1),1);2*ones(size(ct_right.trial,1),1)];
ct_data.trialinfo = clabel;

%% Step 2: Extract warping sources, which contain warping signals
cfg              = [];
cfg.method       = 'runica';
cfg.runica.pca   = 30;               % Let's get 30 to save time
warpsources      = ft_componentanalysis(cfg ,ct_data);


%% Step 3: Time frequency analysis of warping components
% Toolbox configuration
cfg              = [];               % toolbox configuration structure
cfg.time         = [1 2];            % time window of interest (1 to 2s)
cfg.warpfreqs    = [10 10];          % frequency range of interest for brain time (10 to 10 Hz; usually a band)
cfg.correct1f    = 'yes';            % apply 1/f correction after FFT, for plotting purposes only
cfg.nwarpsources  = 10;              % consider only the 10 best warping sources
cfg.rankmethod   = 'maxpow';         % sort warping sources by power in frequency range of interest.
                                     % Alternative: 'templatetopo' (see tutorial 4)
cfg.cutmethod    = 'consistenttime'; % 'cutartefact' or 'consistenttime' See "help bt_analyzecarriers"
                                     % or our paper for details

% Fieldtrip configuration
cfgFT            = [];               % FieldTrip configuration structure, input to ft_freqanalysis
cfgFT.method     = 'wavelet';        % frequency analysis method (see help ft_freqanalysis for options)
cfgFT.width      = 5;                % number of wavelet cycles
cfgFT.foi        = 2:30;             % frequency range for FFT
% cfgFT.time       = cfg.time        % the time variable is automatically fetched from cfg.time

[fft_source_constime]   = bt_analyzesources(cfg,cfgFT,warpsources);

cfg.cutmethod   = 'cutartefact';

[fft_source_cutarte]   = bt_analyzesources(cfg,cfgFT,warpsources);

load layout_tutorial
cfg                      = [];
cfg.layout               = layout;           % load template for topography plotting (not always needed)
[bt_source_constime]     = bt_selectsource(cfg,fft_source_constime,warpsources);
[bt_source_cutarte]      = bt_selectsource(cfg,fft_source_cutarte,warpsources);


%% Check data-to-cycle allocation
cfg.warpmethod   = 'stationary';     % 'stationary': the toolbox warps the data using a stationary
                                     % sinusoid at the warping frequency (default).
                                     % 'waveshape': the toolbox warps to the average waveshape for
                                     % the warping frequency.

cfg.phasemethod  = 'FFT';            % 'FFT': the phase dynamics of the warping signal
                                     % as estimated by ft_freqanalysis (default).
                                     % 'GED': the phase dynamics of the warping signal 
                                     % as estimated by Generalized Eigendecomposition

bt_checkallocation(cfg, ct_data, bt_source_constime, bt_source_cutarte);