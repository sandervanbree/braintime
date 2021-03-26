%%% In tutorial 5 we will make a template topography that can be 
%%% used to bias the ranking of warping sources.
%%% In addition, we will take this chance to look at another way of
%%% analyzing TGMs: 'cutartefact'.

%% Generate a template topography
% Create a template topography, which will be saved in the toolbox
% topography folder. Please draw two boxes with your cursor, one over each
% parietal cortex. This is where we've simulated and therefore predict the
% optimal warping signal to reside.

load layout_tutorial
cfg.layout = layout;
bt_templatetopo(cfg);

%% This section is unchanged from tutorial 1
load dipolesim_tutorial;

% The Brain Time Toolbox requires FieldTrip formatted datastructures, which
% contain both classes of data, labelled 1 and 2 in the trialinfo field
cfg               = [];
ct_data           = ft_appenddata(cfg,ct_left,ct_right);
clabel            = [ones(size(ct_left.trial,1),1);2*ones(size(ct_right.trial,1),1)];
ct_data.trialinfo = clabel;

% Independent component analysis
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

% Fieldtrip configuration
cfgFT            = [];               
cfgFT.method     = 'wavelet';        
cfgFT.width      = 5;                
cfgFT.foi        = 2:30;             

[fft_sources]   = bt_analyzesources(cfg,cfgFT,warpsources);

%% Select a warping source
% Select a warping source with high power at 10 Hz, and parietal sources
load layout_tutorial
cfg              = [];
cfg.layout       = layout;           % load template for topography plotting (not always needed)
[bt_source]      = bt_selectsource(cfg,fft_sources,warpsources);

%% clock time to brain time
cfg              = [];
cfg.bt_srate     = 200;              % specify the sampling rate of bt_data
cfg.removecomp   = 'no';             % remove component when using brain time warped data outside the toolbox to avoid circularity
[bt_warpeddata]  = bt_clocktobrain(cfg,ct_data,bt_source);

% cut ct_data to the same time window
cfg        = [];
cfg.toilim = [bt_warpeddata.toi(1) bt_warpeddata.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

% Feel free to use this output in tutorial 3
save tutorial 5