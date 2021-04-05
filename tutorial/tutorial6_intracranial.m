%%% In tutorial 6 we will brain time warp intracranial data. The data
%%% consists of spike series and LFP obtained from the rat hippocampus.
%%% ASK SIMON FOR SHORT METHODOLOGICAL DETAILS
%%%
%%% There are no limitations to the type of electrophysiological data that 
%%% can be used in the toolbox. All that is required is a FieldTrip
%%% formatted datastructure with the to-be-warped data (clock time data),
%%% and a FieldTrip structure with warping sources. Thus, EEG, MEG, and
%%% intracranial data all work.

% Load two classes of data and carrier LFP (see tutorial folder)
load c1_data_tut5 % contains spike and LFP data for c1
load c2_data_tut5 % contains spike and LFP data for c2

% Label the classes (.trialinfo) and combine spike data
c1_spikes.trialinfo   = ones(size(c1_spikes.trial,1),1);
c2_spikes.trialinfo   = 2*ones(size(c2_spikes.trial,1),1);
cfg                   = [];
ct_data               = ft_appenddata(cfg,c1_spikes,c2_spikes);

% Label the classes (.trialinfo) and combine LFP data
c1_carrier.trialinfo  = ones(size(c1_carrier.trial,1),1);
c2_carrier.trialinfo  = 2*ones(size(c2_carrier.trial,1),1);
cfg                   = [];
carrier_data          = ft_appenddata(cfg,c1_carrier,c2_carrier);

%% Perform FFT over channels to enable sorting by power and to enable phase extraction
cfg = [];
cfg.time         = [0 1.5];          % time window of interest
cfg.fft          = [2 30];           % frequency range for the FFT
cfg.foi          = [4 8];            % for this dataset we will look at the theta range
cfg.waveletwidth = 5;                % wavelet width in number of cycles
cfg.Ntopchans     = 10;               % consider only the 10 best components
cfg.cutmethod    = 'cutartefact'; % 'cutartefact' or 'consistenttime' See "help bt_analyzecarriers" or our paper for details
cfg.sortmethod   = 'maxpow';         % sort by power in frequency range of interest. Alternative: 'templatetopo' (see tutorial 4)
[fft_channels]    = bt_analyzechannels(cfg,carrier_data);

%% Choose a carrier
% Choose the first component's carrier (8 Hz)
cfg              = [];
[bt_carrier]     = bt_choosecarrier(cfg,fft_channels,carrier_data);

%% Warp original clock time data to brain time
cfg              = [];
[bt_struc]       = bt_clocktobrain(cfg,ct_data,bt_carrier);

% cut ct_data to the same window
cfg        = [];
cfg.toilim = [bt_struc.toi(1) bt_struc.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

%% Save results for tutorial 5
save tutorial5_output bt_struc ct_data 

% Feel free to enter these data into tutorial 2 to test for recurrence with this intracranial data.