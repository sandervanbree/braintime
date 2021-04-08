%%% In tutorial 5 we will change several parameters in the brain time
%%% pipeline. For one, we will make a template topography that can be 
%%% used to change the ranking of warping sources to activity in regions of
%%% interest. Another change is that we apply the cutartefact method
%%% during bt_analyzesources, to cope with an artefact that may alter
%%% subsequent analyses. We will also change the warping and phase method
%%% used in bt_clocktobrain. Several tips will be included.

%% Generate a template topography
% Create a template topography, which will be saved in the toolbox
% topography folder. Please draw two boxes with your cursor, one over each
% parietal cortex. To continue, click in one of the drawn boxes.
% Parietal regions is where we've simulated and therefore predict the
% optimal warping signal to reside.
% Template creation is only relevant when using independent component
% analysis (ICA) components as your potential warping sources.

load layout_tutorial
cfg.layout = layout;
bt_templatetopo(cfg);

%% This section is unchanged from tutorial 1
load dipolesim_tutorial;

cfg               = [];
ct_data           = ft_appenddata(cfg,ct_left,ct_right);
clabel            = [ones(size(ct_left.trial,1),1);2*ones(size(ct_right.trial,1),1)];
ct_data.trialinfo = clabel;  % Be sure to overwrite your data's default
                             % trialinfo structure with a vector of 
                             % ones and twos, representing your two classes
                             % (conditions of interest).

cfg              = [];
cfg.method       = 'runica';
cfg.runica.pca   = 30;               % Let's get 30 to save time
warpsources      = ft_componentanalysis(cfg ,ct_data);

% Step 2: Extract warping sources, which contain warping signals
cfg              = [];              
cfg.time         = [1 2];            
cfg.warpfreqs    = [10 10];          
cfg.correct1f    = 'yes';          
cfg.nwarpsources  = 10;           
cfg.rankmethod   = 'templatetopo';  % With 'templatetopo', the toolbox
                                    % biases the warping source sorting to
                                    % the topography drawn using
                                    % bt_templatetopo. The toolbox
                                    % automatically searches for an
                                    % existing file, so be sure to create a
                                    % new template topography when
                                    % performing a new analysis.
                                    %
cfg.cutmethod    = 'cutartefact';   % The toolbox warps 0.5s additional
                                    % data before and after the specified
                                    % window of interest that is later
                                    % removed. For the relative trade-offs
                                    % between 'consistenttime' and 
                                    % 'cutartefact', please refer to
                                    % bt_analyzesource as well as the brain
                                    % time toolbox paper.

% Other FieldTrip methods for frequency estimation are available:
% https://www.fieldtriptoolbox.org/reference/ft_freqanalysis/
cfgFT            = [];                 
cfgFT.method     = 'mtmconvol';        % Let's try a method where we convolve
                                       % the data with a complex wavelet.
cfgFT.taper      = 'hanning';          % How to taper the data?
cfgFT.foi        = 2:30;            
cfgFT.t_ftimwin  = 5./cfgFT.foi;       % Length of sliding window in seconds
cfgFT.tapsmofrq  = 2 +(0*(cfgFT.foi)); % Width of frequency smooting in Hz

[fft_sources]   = bt_analyzesources(cfg,cfgFT,warpsources);

%% Select a warping source
% Select a warping source with high power at 10 Hz, and parietal sources
% The ordering of sources is now changed by their correlation to your
% template topography.
load layout_tutorial
cfg              = [];
cfg.layout       = layout;           % load template for topography plotting (not always needed)
[bt_source]      = bt_selectsource(cfg,fft_sources,warpsources);

%% clock time to brain time
cfg              = [];
cfg.removecomp   = 'no';             % remove component when using brain time warped data outside the toolbox to avoid circularity

cfg.warpmethod   = 'waveshape';      % Let's try warping using the average waveshape in the data

cfg.phasemethod  = 'GED';            % Generalized Eigendecomposition can be used to get a 
                                     % holistic phase estimation of the warping frequency
                                     % across the warping sources
                                     
cfg.visualcheck = 'yes';             % With warpmethod 'waveshape', you see an additional step:
                                     % the construction of the template phase signal based
                                     % on the average waveshape (its concatenation, and phase
                                     % estimation).

[bt_warpeddata]  = bt_clocktobrain(cfg,ct_data,bt_source);

% cut ct_data to the same time window
cfg        = [];
cfg.toilim = [bt_warpeddata.toi(1) bt_warpeddata.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

% Feel free to use this output in tutorial 3
save tutorial5