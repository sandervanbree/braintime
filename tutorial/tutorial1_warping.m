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

% Filter the data
cfg = [];
cfg.bpfilter     = 'yes';
cfg.bpfreq       = [2 30];           % Filter between x and y Hz
ct_data          = ft_preprocessing(cfg,ct_data);

%% Run ICA to extract components, one of which will contain our carrier oscillation
cfg              = [];
cfg.method       = 'runica';
cfg.runica.pca   = 30;               % Optional: obtain N component to reduce time
channels         = ft_componentanalysis(cfg ,ct_data);

%% Perform FFT over channels (components) to enable sorting by time frequency characteristics of interest
cfg = [];
cfg.time         = [0 1];            % time window of interest
cfg.fft          = [2 30];           % frequency range for the FFT
cfg.foi          = [6 10];           % frequency range of interest for brain time
cfg.waveletwidth = 5;                % wavelet width in number of cycles
cfg.Ntop         = 10;               % consider only the 10 best components
cfg.cutmethod    = 'consistenttime'; % 'cutartefact' or 'consistenttime' See "help bt_analyzecarriers" or our paper for details
cfg.sortmethod   = 'maxpow';         % sort by power in frequency range of interest. Alternative: 'templatetopo' (see tutorial 4)
[fft_channels]    = bt_analyzechannels(cfg,channels);

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