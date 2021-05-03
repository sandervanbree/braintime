%%% In tutorial 6 we will brain time warp intracranial data. The data
%%% consists of spike series and LFP obtained from the rat hippocampus.
%%%
%%% There are no limitations to the type of electrophysiological data that 
%%% can be used in the toolbox. All that is required is a FieldTrip
%%% formatted datastructure with the to-be-warped data (clock time data),
%%% and a FieldTrip structure with warping sources. Thus, EEG, MEG, and
%%% intracranial data all work.

% Load an example rat's classes of data and carrier LFP (see tutorial folder)
load gridcell_tutorial.mat rat2

c1_data = rat2.c1_data;
c2_data = rat2.c2_data;

% We have local field potential channels and cell spike channels
disp(c1_data.label)
disp(c2_data.label)

% Combine classes and add condition labels as trialinfo
cfg               = [];
ct_data           = ft_appenddata(cfg,c1_data,c2_data);
clabel            = [ones(size(c1_data.trial,2),1);2*ones(size(c2_data.trial,2),1)];
ct_data.trialinfo = clabel;

% The first two channels (LFP data) serve as the warping sources that
% contain theta
cfg = [];
cfg.channel = 1:2;
warpsources = ft_selectdata(cfg,ct_data);

% The last two channels (cell spikes) serve as the clock time data which
% will be warped
cfg.channel = 3:4;
ct_data = ft_selectdata(cfg,ct_data);

% Toolbox configuration
cfg              = [];               
cfg.time         = [-1.25 1.25];     
cfg.warpfreqs    = [4 8];                % Theta oscillations are tightly linked to grid cell functioning
cfg.correct1f    = 'no';           
cfg.nwarpsources  = 4;            
cfg.rankmethod   = 'maxpow';     
cfg.cutmethod    = 'consistenttime';

% Fieldtrip configuration
cfgFT            = [];               
cfgFT.method     = 'wavelet';        
cfgFT.width      = 4;                
cfgFT.foi        = 2:30;         

[fft_sources]   = bt_analyzesources(cfg,cfgFT,warpsources);

% Select source
cfg              = []; % No layout needed as we are working with intracranial data
[bt_source]      = bt_selectsource(cfg,fft_sources,warpsources);

% Clock to brain
cfg              = [];
cfg.warpmethod   = 'waveshape';    % Rodent theta can be quite asymmetric,
                                   % so warping using the average waveshape
                                   % rather than a stationary sinusoid is a
                                   % natural choice.
cfg.phasemethod  = 'fft';          
cfg.visualcheck  = 'on';
cfg.removecomp   = 'no';
[bt_warpeddata]  = bt_clocktobrain(cfg,ct_data,bt_source);

% cut ct_data to the same time window
cfg        = [];
cfg.toilim = [bt_warpeddata.toi(1) bt_warpeddata.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

%% Save results for tutorial 5
save tutorial6_output

% Feel free to enter these data into tutorial 2 to test for recurrence with this intracranial data.