%%% In tutorial 8 we will test whether brain time warping can uncover
%%% nested oscillations. Specifically, we simulate a cross-frequency
%%% coupled model using bt_dipsim_cfc where high frequency oscillators
%%% are amplitude modulated by the frequency oscillators (see  Jensen &
%%% Colgin; 2007) for details). Does brain time warping to one oscillator
%%% reveal a peak at the frequency of the nested oscillator?
%%% In the supplementary material we show our findings on this question,
%%% and this tutorial shows how to test this question yourself.
%%%
%%% This tutorial merges tutorial 1 and 2, and removes previous comments to
%%% keep it brief. Nothing fundamentally changes from the previous
%%% tutorials, we just look for a different result.

%% Step 1: Load and structure data
load dipolesim_cfc;  % Load to-be-warped data with two cross frequency coupling conditions

% Sort out data structure
cfg               = [];
ct_data           = ft_appenddata(cfg,ct_left,ct_right);
clabel            = [ones(size(ct_left.trial,1),1);2*ones(size(ct_right.trial,1),1)];
ct_data.trialinfo = clabel;

%% Step 2: Extract warping sources, which contain warping signals
cfg              = [];
cfg.method       = 'runica';
cfg.runica.pca   = 30;               % Let's get 30 to save time
warpsources      = ft_componentanalysis(cfg ,ct_data);

%% Step 3: Time frequency analysis of warping sources
% Toolbox configuration
cfg              = [];               
cfg.time         = [1 2];        % Let us slightly increase the time window for 
                                     % improved spectral resolution (0.5 to 2.5s)
cfg.warpfreqs    = [10 10];          % You may warp to either the low frequency
% cfg.warpfreqs    = [17 17];        % Or the nested high frequency
cfg.correct1f    = 'yes';            
cfg.nwarpsources  = 10;              
cfg.rankmethod   = 'maxpow';         
                                     
cfg.cutmethod    = 'consistenttime'; 
                                     
% Fieldtrip configuration
cfgFT            = [];               
cfgFT.method     = 'wavelet';        
cfgFT.width      = 5;                
cfgFT.foi        = 2:30;             
% cfgFT.time       = cfg.time        

[fft_sources]   = bt_analyzesources(cfg,cfgFT,warpsources);

%% Step 4: Select a warping source
% Select a warping source with high power at the high or low simulated frequency
load layout_tutorial
cfg              = [];
cfg.layout       = layout;   
cfg.quickselect  = 1;
[bt_source]      = bt_selectsource(cfg,fft_sources,warpsources);

%% Step 5: Warp clock time to brain time
cfg              = [];
cfg.removecomp   = 'yes';            % Let us remove the component this time
cfg.warpmethod   = 'sinusoid';       
cfg.phasemethod  = 'fft';
cfg.visualcheck  = 'on';             % It is always good to check whether warping looks reasonable
[bt_warpeddata]  = bt_clocktobrain(cfg,ct_data,bt_source);

% cut ct_data to the same time window
cfg        = [];
cfg.toilim = [bt_warpeddata.toi(1) bt_warpeddata.toi(2)];
ct_data    = ft_redefinetrial(cfg, ct_data);

% Save results (OPTIONAL; we continue here)
% save tutorial8_output bt_warpeddata ct_data warpsources

% Your data is now in brain time. Now let us continue with tutorial 2 to
% figure out whether the nested oscillation shows up.

%% Tutorial 2
% Timelock the data
cfg = [];
cfg.keeptrials     = 'yes';
cfg.removemean     = 'no';
ct_data            = ft_timelockanalysis(cfg,ct_data);
bt_data            = ft_timelockanalysis(cfg,bt_warpeddata.data);

% Extract classification labels (clabels)
clabel             = bt_warpeddata.clabel;

%% Create Clock and Brain Time TGMs
% Use MVPA Light to generate time generalization matrices (TGM)
cfg_mv.classifier  = 'lda';     
cfg_mv.metric      = 'acc';     
cfg_mv.repeat      = 3;         
cfg_mv.cv          = 'kfold'; 
cfg_mv.k           = 5;         
[~, ct_mv]         = mv_classify_timextime(cfg_mv, ct_data.trial, clabel);
[~, bt_mv]         = mv_classify_timextime(cfg_mv, bt_data.trial, clabel); 

%% Quantify TGM periodicity (compare clock and brain time)
cfg = [];
cfg.bt_warpeddata   = bt_warpeddata;              
cfg.figure          = 'yes';
cfg.method          = 'tgm';                      % For Cross-frequency coupling, TGM is the desired method
                                                  % because AC tends to drown out all but the
                                                  % most dominant frequency
cfg.periodicityfoi   = [1 30];                    % Make sure to include both frequencies in this range
cfg.refdimension    = 'clocktime';                
ct_quant         = bt_quantify(cfg,ct_mv);        
title('Clock time periodicity');

% To make the comparison fair visually, we need to equalize the y-axes for clock and brain time
yl = get(gca,'YLim');

cfg.refdimension    = 'braintime';                
bt_quant         = bt_quantify(cfg,bt_mv);        
title('Brain time periodicity');
ylim(yl); % Adapy y-axis

% Compare the power at the nested frequency between clock and brain time.
% You could also continue with tutorial 3 to include statistics