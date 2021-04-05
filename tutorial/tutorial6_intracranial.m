%% NOTE TO DEVS: This tutorial requires some patching up. Note to self:
%%% Ask simon to repeat methodological details of the rodent data
%%%
%%% In tutorial 6 we will brain time warp intracranial data. The data
%%% consists of spike series and LFP obtained from the rat hippocampus.
%%%
%%% There are no limitations to the type of electrophysiological data that 
%%% can be used in the toolbox. All that is required is a FieldTrip
%%% formatted datastructure with the to-be-warped data (clock time data),
%%% and a FieldTrip structure with warping sources. Thus, EEG, MEG, and
%%% intracranial data all work.

% Load two classes of data and carrier LFP (see tutorial folder)
load gridcell_tutorial.mat

% NOTE TO DEVS: Include information about channels
disp(c1_data.label)
disp(c2_data.label)

% Combine classes and add condition labels as trialinfo
cfg               = [];
ct_data           = ft_appenddata(cfg,c1_data,c2_data);
clabel            = [ones(size(c1_data.trial,2),1);2*ones(size(c2_data.trial,2),1)];
ct_data.trialinfo = clabel;

% The clock time data also serves as the warping source
warpsources      = ct_data;

% But only the first two channels
cfg = [];
cfg.channel = 1:2;
warpsources = ft_selectdata(cfg,warpsources);

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

% No layout needed as we are working with intracranial data
cfg              = [];
cfg.warpmethod   = 'waveshape';    % Rodent theta can be quite asymmetric,
                                   % so warping using the average waveshape
                                   % rather than a stationary sinusoid is a
                                   % natural choice.
cfg.phasemethod  = 'FFT';          
cfg.visualcheck  = 'on';
cfg.removecomp   = 'no';
[bt_warpeddata]  = bt_clocktobrain(cfg,ct_data,bt_source);

% cut ct_data to the same window
cfg        = [];
cfg.toilim = [bt_struc.toi(1) bt_struc.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

%% Save results for tutorial 5
save tutorial6_output

% Feel free to enter these data into tutorial 2 to test for recurrence with this intracranial data.