%%% In tutorial 6 we will apply an alternative way of finding the carrier
%%% oscillation. As described by Cohen, 2017; DOI: 10.7554/eLife.21792,
%%% Generalized eigendecomposition (GED) can optimize carrier selection
%%% by finding an optimal separation between frequencies of interest
%%% and frequences of no interest. For GED, the analyzing of the data
%%% and choosing of a carrier oscillation are combined into one function,
%%% 'bt_GEDanalyzechoose'. This output can be used in 'bt_clocktobrain'.
%%%
%%% TO DEVS: bt_GEDanalyzechoose is currently not working yet as intended.
%%% Help needed. See the function for more details.

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

% Preprocessed ct_data can serve as channels
channels = ct_data;

%% Generalized eigendecomposition as a method to analyze channels and choose a carrier oscillation
cfg = [];
cfg.time         = [0 1.2];          % time window of interest
cfg.foi          = [6 10];           % frequency range of interest for brain time
cfg.fwhm         = 5;                % full width half maximum for GED
cfg.Ntopchan     = 10;               % consider only the 10 best components
cfg.removecomp   = 'no';             % specify component removal here
cfg.cutmethod    = 'consistenttime'; % cut method
[bt_carrier] = bt_GEDanalyzechoose(cfg,channels); % when cfg.removecomp = 'yes', specify an additional "channels" output to receive the channel data with component removed

%% Warp original clock time data to brain time
cfg              = [];
cfg.btsrate      = 128;              % determine sampling rate of bt data
[bt_struc]       = bt_clocktobrain(cfg,ct_data,bt_carrier);

% cut ct_data to the same window
cfg        = [];
cfg.toilim = [bt_struc.toi(1) bt_struc.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

%% Save results for tutorial 6
save tutorial6_output bt_struc ct_data 

% Feel free to enter these data into tutorial 2 to test for recurrence with this data.