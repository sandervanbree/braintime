%%% In tutorial 1 we will transform data from clock time (input) to
%%% brain time (output). The tutorial loads two classes of 3 second
%%% simulated data (-1s to 2s) with a 8 Hz recurring pattern.
%%% We will extract and warp the first second (0 to 1s) of data.

%%% After combining both classes of data and applying preprocessing,
%%% we will perform ICA on the combined data to find a component with
%%% our "carrier" neural oscillation. The clock time data will be warped
%%% to the phase of this carrier, transforming clock to brain time.
%%% ct = clock time, bt = brain time. See Github for a full glossary.

% Load two classes of data (see tutorial folder)
load c1_data
load c2_data

% Combine classes and add condition labels
cfg               = [];
ct_data           = ft_appenddata(cfg,c1_data,c2_data);
clabel            = [ones(size(c1_data.trial,1),1);2*ones(size(c2_data.trial,1),1)];
ct_data.trialinfo = clabel;

%% Run ICA to extract components, one of which will contain our carrier oscillation
cfg              = [];
cfg.method       = 'runica';
cfg.runica.pca   = 30;               % Optional: obtain N component to reduce time
channels         = ft_componentanalysis(cfg ,ct_data);

%% Perform FFT over channels (components) to enable sorting by time frequency characteristics of interest

% Toolbox configuration
cfg              = [];               % Toolbox configuration structure
cfg.time         = [0 1];            % time window of interest
cfg.warpfreqs    = [6 10];           % frequency range of interest for brain time
cfg.correct1f    = 'yes';            % apply 1/f correction after FFT, for plotting purposes only
cfg.ntopchans    = 10;               % consider only the 10 best components
cfg.sortmethod   = 'maxpow';         % sort by power in frequency range of interest. Alternative: 'templatetopo' (see tutorial 4)
cfg.cutmethod    = 'cutartefact';    % 'cutartefact' or 'consistenttime' See "help bt_analyzecarriers" or our paper for details

% Fieldtrip configuration
cfgFT            = [];               % FieldTrip configuration structure, input to ft_freqanalysis
cfgFT.method     = 'wavelet';        % Frequency analysis method (see ft_freqanalyis)
cfgFT.width      = 5;                % Number of wavelet cycles
cfgFT.foi        = 2:30;             % Frequency range for FFT
% cfgFT.time       = cfg.time        % The time variable is automatically grabbed from cfg.time

[fft_channels]   = bt_analyzechannels(cfg,cfgFT,channels);

%% Choose a carrier
% Choose the first component's carrier (8 Hz)
load layout_tutorial
cfg              = [];
cfg.layout       = layout;           % load template for topography plotting
[bt_carrier]     = bt_choosecarrier(cfg,fft_channels,channels);

%% Warp original clock time data to brain time
cfg              = [];
cfg.btsrate      = 128;              % determine sampling rate of bt data
cfg.removecomp   = 'no';             % remove component when using brain time warped data outside the toolbox to avoid circularity
[bt_struc]       = bt_clocktobrain(cfg,ct_data,bt_carrier);

% cut ct_data to the same window
cfg        = [];
cfg.toilim = [bt_struc.toi(1) bt_struc.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

%% Save results for tutorial 2
save tutorial1_output bt_struc ct_data