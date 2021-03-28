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
% to facilitate time frequency analyses
load dipolesim_tutorial;

% The Brain Time Toolbox requires FieldTrip formatted datastructures, which
% contain both classes of data, labelled 1 and 2 in the trialinfo field
cfg               = [];
ct_data           = ft_appenddata(cfg,ct_left,ct_right);
clabel            = [ones(size(ct_left.trial,1),1);2*ones(size(ct_right.trial,1),1)];
ct_data.trialinfo = clabel;

% It is important that the data is not or minimally high pass
% filtered, as this can introduce strong artifacts in later steps in the
% toolbox. See van Driel et al., 2021; https://doi.org/10.1016/j.jneumeth.2021.109080
% for details and recommendations.

%% Step 2: Extract warping sources, which contain warping signals
% "Warping source": The data that contains your warping signal. These can
% be ICA components, source localized data, LFP from macrowires, etc.
% "Warping signal": A signal of interest within warping sources. The
% clock time data will be warped according to this signal's dynamics.
cfg              = [];
cfg.method       = 'runica';
cfg.runica.pca   = 30;               % Let's get 30 to save time
warpsources      = ft_componentanalysis(cfg ,ct_data);

% It is recommended to save your warpsources (ICA output), to replicate
% your results at a later stage:
% save warpsources

%% Step 3: Time frequency analysis of warping components
% Toolbox configuration
cfg              = [];               % toolbox configuration structure
cfg.time         = [1 2];            % time window of interest (1 to 2s)
cfg.warpfreqs    = [10 10];          % frequency range of interest for brain time (10 to 10 Hz; usually a band)
cfg.correct1f    = 'yes';            % apply 1/f correction after FFT, for plotting purposes only
cfg.nwarpsources  = 10;              % consider only the 10 best warping sources
cfg.rankmethod   = 'maxpow';         % sort warping sources by power in frequency range of interest.
                                     %  Alternative: 'templatetopo' (see tutorial 4)
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
cfg.bt_srate     = ct_data.fsample;  % specify the sampling rate of bt_data
cfg.removecomp   = 'no';             % remove component when using brain time warped data outside
                                     % the toolbox to avoid circularity
[bt_warpeddata]  = bt_clocktobrain(cfg,ct_data,bt_source);

% cut ct_data to the same time window
cfg        = [];
cfg.toilim = [bt_warpeddata.toi(1) bt_warpeddata.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

% Save results for tutorial 2
save tutorial1_output bt_warpeddata ct_data

% Your data is now in brain time. This means the original data is rescaled
% based on the phase dynamics of your warping signal. You can now use your
% transformed data in your own analyses (though ensure cfg.removecomp was
% on 'yes' to avoid circularity). You may also test whether brain time
% warping has unveiled periodic patterns in your data by applying
% classification. To see how this is done, check out the next tutorials.