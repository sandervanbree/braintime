%%% In tutorial 4 we will make a template topography that can be 
%%% used for choosing the most suitable brain time ICA component.
%%% We will also try different parameters in the pipeline

% Create template topography (saved in topography folder)
% Draw the box over frontal areas, where the simulated patterns were induced
load layout_tutorial
cfg.layout = layout;
bt_templatetopo(cfg);

%% This section is unchanged from tutorial 1
% Load two classes of data (see tutorial folder)
load c1_data
load c2_data

% Label the classes (.trialinfo) and combine them
c1_data.trialinfo = ones(size(c1_data.trial,1),1);
c2_data.trialinfo = 2*ones(size(c2_data.trial,1),1);
cfg               = [];
ct_data           = ft_appenddata(cfg,c1_data,c2_data);

% Filter the data
cfg = [];
cfg.bpfilter     = 'yes';
cfg.bpfreq       = [2 30];           % Filter between x and y Hz
ct_data             = ft_preprocessing(cfg,ct_data);

% Run ICA to extract components, one of which will contain the designated brain time oscillation
cfg              = [];
cfg.method       = 'runica';
cfg.runica.pca   = 30;               % Optional: obtain N component, to reduce time
comp             = ft_componentanalysis(cfg ,ct_data);

%% Perform FFT over components to enable sorting by power and to enable phase extraction
% This time, we will choose the cutmethod 'cutartefact', which gets rid of 
% the artefact at the first cycle, but causes different trial durations,
% resulting in reduced recurrence.

% In addition, we will choose sortmethod 'temptopo', which sorts components
% by their ranking on two factors: power and correlation to your template
% topography.

cfg = [];
cfg.time         = [0 1];            % time window of interest
cfg.fft          = [2 30];           % frequency range for the FFT
cfg.foi          = [6 10];           % frequency range of interest for brain time
cfg.waveletwidth = 5;                % wavelet width in number of cycles
cfg.topcomp      = 10;               % number of components to be considered
cfg.cutmethod    = 'cutartefact';    % 'cutartefact' or 'consistenttime' See "help bt_analyzecomps" or our paper for details
cfg.sortmethod   = 'templatetopo';   % sort by power in frequency range of interest. Alternative: 'templatetopo' (see tutorial 2)
[fft_comp]       = bt_analyzecomps(cfg,comp);

%% Designate component frequency as brain time
% Choose component with a frontal topography at 8 Hz
load layout_tutorial
cfg              = [];
cfg.layout       = layout;           % load template for topography plotting
[bt_comp]        = bt_choosecomp(cfg,fft_comp,comp);

%% Warp original clock time data to brain time
cfg              = [];
cfg.btsrate      = 128;              % determine sampling rate of bt data
cfg.removecomp   = 'yes';            % remove component from data to avoid circularity (see paper)
[bt_struc]        = bt_clocktobrain(cfg,ct_data,bt_comp);

% cut ct_data to the same window
cfg        = [];
cfg.toilim = [bt_struc.toi(1) bt_struc.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

%% Save results
save tutorial4_output bt_struc ct_data

% Feel free to enter these data into tutorial 2 to test for recurrence with the new parameters.