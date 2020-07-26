%%% In tutorial 5 we will brain time warp intracranial data, based on
%%% recordings with a Behnke Fried electrode. In this tutorial, the clock
%%% time data and carrier oscillation come from different recording
%%% sources. Clock time data are spike trains convolved with a smoothing
%%% kernel. The channel data with electable carrier oscillations are LFP
%%% recordings. This differs from tutorial 1, where clock time data
%%% and carrier channel data came from the same data structure.
%%%
%%% TO DEVS: I have found intracranial data that is more suited (i.e. it
%%% shows recurrence only in brain time), but it is too large for Github.
%%% As soon as I found a way to properly downsample, I will replace
%%% c1_data_tut5 and c2_data_tut5.

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

% Filter the data
cfg = [];
cfg.bpfilter          = 'yes';
cfg.bpfreq            = [2 30]; % Filter between x and y Hz
ct_data               = ft_preprocessing(cfg,ct_data);
carrier_data          = ft_preprocessing(cfg,carrier_data);

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
cfg.btsrate      = 200;              % determine sampling rate of bt data
[bt_struc]       = bt_clocktobrain(cfg,ct_data,bt_carrier);

% cut ct_data to the same window
cfg        = [];
cfg.toilim = [bt_struc.toi(1) bt_struc.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

%% Save results for tutorial 5
save tutorial5_output bt_struc ct_data 

% Feel free to enter these data into tutorial 2 to test for recurrence with this intracranial data.