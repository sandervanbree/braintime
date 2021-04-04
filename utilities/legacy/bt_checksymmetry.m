function [symm] = bt_checksymmetry(config,warpsources)
% This is a legacy function that is no longer part of the Brain Time
% Toolbox, but may come in useful. This function test the asymmetry of all
% warping sources. It extracts the average waveshape and asymmetry indices
% for each desired frequency. It returns a ranking of warping sources
% from most to least symmetric, with t-statistics and waveshapes.
%
% Use as:
%
% [symm] = bt_checksymmetry(cfg, warpsources)
%
% Input Arguments:
% config             %
%                    %
% - foi              % Range of integer frequencies over which asymmetry is
%                    % tested for each warping source
%                    %
% Output:            %
% symm               % Output structure with asymmetry indices and
%                    % waveshapes for each frequency (rows) and each
%                    % warping source (columns; ranked from most to least
%                    % symmetric.

%% Get basic info
foi            = config.foi;                                 % Frequencies for which asymmetry is tested
ncycles        = 2;                                          % Hard-code number of cycles to display (visualization only)
wsources       = warpsources;                                % Rename for brevity
nwsources      = size(wsources.trial{1},1);                  % Fetch number of warping sources

%% Data handling
% Slice data to time window of interest
if isfield(config,'toi') == 1
    cfg             = [];
    cfg.latency     = config.toi;
    wsources        = ft_selectdata(cfg,wsources);
end

% Get average ASI and waveshape across trials
[asymmidx,asymmidx_t,wavshap] = bt_sourcesymm(wsources,foi,ncycles);

% Compute average across warping sources
wavshap_f = squeeze(mean(wavshap,2));

figure; hold on; bt_figure;
subplot(1,2,1)
bt_rainplot_legacy(asymmidx,[],50,100,foi);      % Plot raincloud plot for each frequency
title('Asymmetry index across warping sources');
subplot(1,2,2)
bt_wavplot_legacy(wavshap_f,ncycles,foi);        % Plot average waveshapes next to it
title('Average waveshape across warping sources');

% Create output structure with all variables, ordered by most to least symmetric
[~,ws_rank]        = sort(abs(asymmidx),2,'ascend');
symm.wsources      = ws_rank;                          % Order warping sources index by most to least symmetric
symm.f             = foi;                              % Tested frequencies
% For every frequency
for f = 1:numel(foi)
    symm.asymmidx(f,:)      = asymmidx(f,ws_rank(f,:));    % Add and order asymmetry index
    symm.asymmidx_t(f,:)    = asymmidx_t(f,ws_rank(f,:));  % Add and order asymmetry index' t-values
    
    % For every timepoint of the waveshape
    for t = 1:size(wavshap,3)
        tmp = wavshap(f,:,t);
        symm.wavshape(f,:,t) = tmp(ws_rank(f,:));              % Add and order waveshape
    end
end

symm.outputformat  = ['Frequencies ',num2str(foi(1)),' to ',num2str(foi(end)),' Hz (rows) and warping source 1 to ',num2str(nwsources),' (columns), sorted by most to least symmetric warping source.'];

disp('The output structure is organized in the following way:');
disp(symm.outputformat);
end
