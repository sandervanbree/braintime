%%% In tutorial 1 we will transform data from clock time (input) to
%%% brain time (output). The tutorial loads two classes of simulated data,
%%% going from -1 to 2 seconds. The data contains an 8-Hz recurring
%%% pattern. We will extract and warp the middle portion of data (0 to 1s).
%%%
%%% After combining both classes of data and applying preprocessing,
%%% we will perform ICA on the combined data to find and rank warping
%%% sources according to their time frequency characteristics of interest.
%%% Next, we will select the optimal warping source. The warping signal
%%% in the selected warping source will be used to transform the clock time
%%% data into brain time.
%%%
%%% ct = clock time, bt = brain time. See Github for a full glossary.

% Load two classes of data (see tutorial folder)
% Ensure that the data
% (1) is longer than your time window of interest to facilitate time
% frequency analysis
% (2) has a trialinfo field that can be used to extract class labels of
% your two conditions of interest
load c1_data
load c2_data

% It is very important that the data is not or minimally high pass
% filtered, as this can introduce strong artifacts in later steps in the
% toolbox (see van Driel et al., 2021; https://doi.org/10.1101/530220)

% Combine classes and add condition labels as trialinfo
cfg               = [];
ct_data           = ft_appenddata(cfg,c1_data,c2_data);
clabel            = [ones(size(c1_data.trial,1),1);2*ones(size(c2_data.trial,1),1)];
ct_data.trialinfo = clabel;

%% Run ICA to extract warping sources, one of which will contain our warping signal
cfg              = [];
cfg.method       = 'runica';
cfg.runica.pca   = 30;               % Let's get 30 to save time
warpsources      = ft_componentanalysis(cfg ,ct_data);

% It is recommended to save your warpsources (ICA output),
% to replicate your results at a later stage.
save warpsources

%% Perform FFT over warping sources (components) to enable sorting by time frequency characteristics of interest

% Toolbox configuration
cfg              = [];               % toolbox configuration structure
cfg.time         = [0 1];            % time window of interest (0 to 1s)
cfg.warpfreqs    = [6 10];           % frequency range of interest for brain time (6 to 10Hz)
cfg.correct1f    = 'yes';            % apply 1/f correction after FFT, for plotting purposes only
cfg.nwarpsources  = 10;              % consider only the 10 best warping sources
cfg.rankmethod   = 'maxpow';         % sort warping sources by power in frequency range of interest. Alternative: 'templatetopo' (see tutorial 4)
cfg.cutmethod    = 'consistenttime';    % 'cutartefact' or 'consistenttime' See "help bt_analyzecarriers" or our paper for details

% Fieldtrip configuration
cfgFT            = [];               % FieldTrip configuration structure, input to ft_freqanalysis
cfgFT.method     = 'wavelet';        % frequency analysis method (see help ft_freqanalysis for options)
cfgFT.width      = 5;                % number of wavelet cycles
cfgFT.foi        = 2:30;             % frequency range for FFT
% cfgFT.time       = cfg.time        % the time variable is automatically grabbed from cfg.time

[fft_sources]   = bt_analyzesources(cfg,cfgFT,warpsources);

%% Select a warping source
% For this tutorial, select the first warping source, which contains a warping signal at 8 Hz
load layout_tutorial
cfg              = [];
cfg.layout       = layout;           % load template for topography plotting
[bt_source]      = bt_selectsource(cfg,fft_sources,warpsources);

%% Warp original clock time data to brain time
cfg              = [];
cfg.btsrate      = 200;              % specify the sampling rate of bt_data
cfg.removecomp   = 'no';             % remove component when using brain time warped data outside the toolbox to avoid circularity
[bt_warpeddata]  = bt_clocktobrain(cfg,ct_data,bt_source);

% cut ct_data to the same window
cfg        = [];
cfg.toilim = [bt_warpeddata.toi(1) bt_warpeddata.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

%% Save results for tutorial 2
save tutorial1_output bt_warpeddata ct_data