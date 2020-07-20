%%% In tutorial 7: data-to-cycle allocation. Not documented yet.

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
cfg.cutmethod    = 'consistenttime'; % cutmethod 1
cfg.sortmethod   = 'maxpow';         % sort by power in frequency range of interest. Alternative: 'templatetopo' (see tutorial 4)
[fft_channels_1]    = bt_analyzechannels(cfg,channels);

cfg.cutmethod    = 'cutartefact';    % cutmethod 2
[fft_channels_2]    = bt_analyzechannels(cfg,channels);

%% Choose a carrier
% Choose the first component's carrier (8 Hz)
load layout_tutorial
cfg              = [];
cfg.layout       = layout;           % load template for topography plotting
[bt_carrier_1]     = bt_choosecarrier(cfg,fft_channels_1,channels);
[bt_carrier_2]     = bt_choosecarrier(cfg,fft_channels_2,channels);

%% Check data-to-cycle allocation
cfg              = [];
cfg.btsrate      = 128;              % determine sampling rate of bt data
bt_checkallocation(cfg, ct_data, bt_carrier_1, bt_carrier_2);