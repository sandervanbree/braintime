%%% In this script we will use the legacy function bt_checksymmetry to
%%% test the asymmetry of all warping sources. The function extracts the
%%% average waveshape and asymmetry indices for each desired frequency.
%%% It returns a ranking of warping sources from most to least symmetric,
%%% with t-statistics and waveshapes.

%% Load some warping sources
load tutorial1_output warpsources

%% Input parameters
cfg         = [];
cfg.foi     = 8:1:12;          % Range of integer frequencies for which to test asymmetry
[symm]      = bt_checksymmetry(cfg, warpsources);

% Importantly, dots in this raincloud represent the average asymmetry per
% warping source. The primary (live) bt_rainplot function's dots represent the
% average asymmetry per trial.

%% Results extraction
disp(symm.outputformat);       % Each output cell is structured according to
                                 % frequencies (rows) and warping source
                                 % indices (columns).

disp(symm.wsources);           % Per frequency (rows) each warping source is ranked (columns)
disp(symm.f);                  % Tested frequencies
disp(symm.asymmidx);           % Asymmetry indices 
disp(symm.asymmidx_t);         % Asymmetry t-statistics
disp(symm.wavshape);           % Average waveshapes

