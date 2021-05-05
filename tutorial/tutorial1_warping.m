%%% In tutorial 1 we will transform data from clock time (input) to
%%% brain time (output). The tutorial loads two classes of simulated data,
%%% generated using dipole simulation. We will brain time warp the middle
%%% segment of the dataset: 1 to 2s.
%%% 
%%% The simulation comprises a simple rhythmic attentional model, where
%%% attention is oriented to the left (class 1) or right (class 2)
%%% hemifield. The dipoles of interest oscillate at 10 Hz, but start
%%% with a different phase between trials, and show frequency drift within
%%% trials. For details on this simulation, please refer to:
%%% (.../dipolesimulation/bt_dipsim.m).
%%%
%%% ct = clock time, bt = brain time. See Github for a full glossary.

%% Step 1: Load and structure data
% Ensure your data is longer than your time window of interest,
% to facilitate time frequency analyses.
load dipolesim_tutorial;  % Load to-be-warped data with two conditions

% The Brain Time Toolbox requires the data to be in a FieldTrip format
cfg               = [];
ct_data           = ft_appenddata(cfg,ct_left,ct_right);

% In case the second operation of the toolbox will be used (periodicity
% analysis; tutorial 2 and 3), please enter the class labels in the data's trialinfo
% field. If only he first operation of the toolbox will be
% used (brain time warping), this is not necessary.
clabel            = [ones(size(ct_left.trial,1),1);2*ones(size(ct_right.trial,1),1)];
ct_data.trialinfo = clabel;

%% Step 2: Extract warping sources, which contain warping signals
% The toolbox requires two electrophysiology data structures: clock time
% data and warping sources. The former is what you want to transform. The 
% latter is what you use to transform. These can be ICA components, source
% localized data, LFP from macrowires, or a channel from the clock time
% data itself.
% Each warping source contains potential warping signals - a frequency of
% interest predicted to carry the studied cognitive process. The
% clock time data will be warped according to this signal's dynamics.

cfg              = [];
cfg.method       = 'runica';
cfg.runica.pca   = 30;               % Let's get 30 to save time
warpsources      = ft_componentanalysis(cfg ,ct_data);

% It is recommended to save your warpsources (ICA output), to replicate
% your results at a later stage:
% save warpsources

%% Step 3: Time frequency analysis of warping sources
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

[fft_sources]   = bt_analyzesources(cfg,cfgFT,warpsources);

%% Step 4: Select a warping source
% Select a warping source with high power at 10 Hz, and parietal sources
load layout_tutorial
cfg              = [];
cfg.layout       = layout;           % load template for topography plotting (not always needed)
[bt_source]      = bt_selectsource(cfg,fft_sources,warpsources);

%% Step 5: Warp clock time to brain time
cfg              = [];
cfg.removecomp   = 'no';             % remove component when using brain time warped data outside
                                     % the toolbox to avoid circularity
cfg.warpmethod   = 'sinusoid';       % 'sinusoid': the toolbox warps the data using a stationary
                                     % sinusoid at the warping frequency (default).
                                     % 'waveshape': the toolbox warps to the average waveshape for
                                     % the warping frequency.

cfg.phasemethod  = 'fft';            % 'fft': the phase dynamics of the warping signal
                                     % as estimated by ft_freqanalysis (default).
                                     % 'ged': the phase dynamics of the warping signal 
                                     % as estimated by Generalized Eigendecomposition
                                     
cfg.visualcheck  = 'on';             % Visualize several steps to check for errors
[bt_warpeddata]  = bt_clocktobrain(cfg,ct_data,bt_source);

% Let's take a look at the output of visualcheck, which illustrates the toolbox's core workings.
% Figure 1: dynamic time warping for two example trials. At the top, you see (1) the unwrapped
% phase of the chosen warping signal (blue; "brain time"), and (2) the unwrapped phase of a
% stationary sinusoid (orange; "clock time"). In the bottom, you see what happens after applying warping: 
% the algorithm attempts to minimize the difference between the two signals.
%
% Figure 2: The main premise of the toolbox is that the warping path between the signals reveals moments
% in each trial when the warping signal - which orchestrates the studied cognitive process - falls out
% of tune with clock time. The toolbox employs the warping path to dynamically resample the data,
% repeating datapoints where brain time slows down relative to clock time, then squeezing
% everything back to the same length. Here you see roughly what the transformation entails on
% a single-trial basis.

% cut ct_data to the same time window
cfg        = [];
cfg.toilim = [bt_warpeddata.toi(1) bt_warpeddata.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

% Save results for tutorial 2
save tutorial1_output bt_warpeddata ct_data warpsources

% Your data is now in brain time. This means the original data is rescaled
% based on the phase dynamics of your warping signal. You can now use your
% transformed data in your own analyses (though ensure cfg.removecomp was
% on 'yes' to avoid circularity). You may also test whether brain time
% warping has unveiled periodic patterns in your data by applying
% classification. To see how this is done, check out the next tutorials.
