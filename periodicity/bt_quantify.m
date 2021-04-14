function [quant] = bt_quantify(config, mv_results)
% Quantify the degree of periodicity in the time generalization matrix
% (TGM), or its diagonal. Here, "periodicity" refers to fluctuating 
% patterns of classification, which reflect fluctuations of the studied
% cognitive process. This function displays the TGM/diag, the
% autocorrelation map (only for TGM). Then, it performs a fast Fourier
% Transform over all the rows and columns (for TGM), or just one row 
% (diagonal). Periodicity quantification is necessary for first level
% statistics (bt_statslevel1).
%
% Use as:
% [bt_quant] = bt_quantify(cfg, mv_results)
% [ct_quant] = bt_quantify(cfg, mv_results) (tip: use name to label output [ct_quant])
%
% Input Arguments:
% config
%   - bt_warpeddata  % brain time data structure as obtained by
%                    % bt_clocktobrain.
%                    %
%   - refdimension   % 'clocktime': find periodicity as a function of
%                    % clock time seconds in the brain time data.
%                    %
%                    % 'braintime': find periodicity as a function of
%                    % cycles of the warped frequency in the brain time
%                    & data.
%                    %
%   - periodicityfoi % '[min max]': Range of periodicity frequencies to
%                    % be quantified and statistically tested.
%                    %
%   - method         % This parameter should only be specified for temporal
%                    % generalization analyis (TGM).
%                    %
%                    % 'tgm': perform periodicity analysis over the tgm
%                    % itself. This method is more sensitive to tgms with
%                    % multiple periodicity frequencies.
%                    %
%                    % 'ac' (default): perform periodicity analysis over the
%                    % autocorrelation map of the tgm. This accentuates
%                    % the primary periodicity frequency in the tgm, but may
%                    % drown out other frequencies.
%                    %
%   - figure         % 'yes' (default) or 'no': display tgm and its ac map.
%                    %
% mv_results         % Results obtained by MVPA Light's mv_classify_timextime
%                    %
% Output:            %
% quant              % Data structure with: 
%                    % MVPA configuration structure, periodicity power
%                    % spectrum, TGM/diag, and method of analysis.

%% Get basic info
toi = config.bt_warpeddata.toi;                       % Start and end time of interest
warpfreq = config.bt_warpeddata.freq;                 % Warped requency (frequency of the carrier)
duration = toi(2)-toi(1);                             % Duration of the time window of interest
mv_cfg = mv_results.cfg;                              % MVPA Light results structure
clabel = config.bt_warpeddata.clabel;                 % Classification labels
mv_perf = mv_results.perf;                            % Temporal generalization matrix
refdimension.dim = config.refdimension;               % Extract reference dimension (clock or brain time)

% Check clabel vector
if sum(unique(clabel))~=3 || isvector(clabel)~= 1
    error(['bt_warpeddata.clabel does not contain a vector of 1''s and 2''s. '...
        'Please overwrite it with a correct classification label vector.']);
end

% Figure out map type
if strcmp(mv_results.function,'mv_classify_across_time')       
    maptype = 'diag';
elseif strcmp(mv_results.function,'mv_classify_timextime')
    maptype = 'tgm';
end

% Lowercase method input
if isfield(config,'method')
    config.method = lower(config.method);
end

% Check whether the method input even exists
if isfield(config,'method')
    if (strcmp(config.method,'ac') || strcmp(config.method,'tgm')) == 0
        error('config.method not recognized. Please change it to ''ac'' or ''tgm''');
    end
end

% For tgm, check whether a method was specified
if strcmp(maptype,'tgm')
    if (isfield(config,'method') && (strcmp(config.method,'tgm') || strcmp(config.method,'ac'))) == 0
    warning_msg1 = (['no cfg.method was submitted, so it is set to ''ac'' by default.'...
        ' See help bt_quantify for details.']);
    config.method = 'ac';
    end
end

% Check compatibility maptype and method
if strcmp(maptype,'diag')
    if isfield(config,'method') && (strcmp(config.method,'tgm') || strcmp(config.method,'ac'))
        error(['Without temporal generalization, cfg.method ''tgm'' and ''ac'' are not applicable.'...
            ' Please remove cfg.method and try again.']);
    else
        config.method = 'diag';
    end
end

% Set up periodicity range over which stats will be applied
if isfield(config,'periodicityfoi')
    powspecrange = config.periodicityfoi(1):config.periodicityfoi(end);
else
    warning_msg2 = (['No periodicity frequency range of interest entered under '...
        'config.periodicityfoi; testing 1 to 30 Hz.']);
    powspecrange = 1:30;
end

if strcmp(refdimension.dim,'braintime') % Test if range is OK
    if powspecrange(1)/warpfreq > 0.5
        answer = questdlg(['The lowest frequency of the specified periodicity range '...
            'is too high to test periodicity at half the warped frequency, which is an'...
            ' important harmonic. Please consider decreasing it before continuing.'], ...
            'Warning', ...
            'Continue','Cancel','Continue');
        
        if strcmp(answer,'Cancel')
            error('Please change cfg.periodicityfoi')
        end
    elseif powspecrange(end)/warpfreq < 2
        
        answer = questdlg(['The highest frequency of the specified recurrence range'...
            ' is too low to test periodicity at double the warped frequency, which is an'...
            ' important harmonic. Please consider increasing it before continuing.'], ...
            'Warning', ...
            'Continue','Cancel','Continue');
        
        if strcmp(answer,'Cancel')
            error('Please change cfg.periodicityfoi')
        end
    end
end

if strcmp(refdimension.dim,'braintime')
    % Adjust powspecrange to be centered on warping frequency
    nrst = nearest(powspecrange,warpfreq);
    diffr = powspecrange(nrst)-warpfreq;
    powspecrange = powspecrange-diffr;
    
    powspecrange = powspecrange/warpfreq; % adjust tested range to warped frequency
    
    refdimension.val = duration*warpfreq; % normalize by cycles in the data
    timevec = config.bt_warpeddata.data.time{1};
elseif strcmp(refdimension.dim,'clocktime')
    refdimension.val = duration; %normalize by seconds in the data
    timevec = linspace(toi(1),toi(2),size(mv_perf,1));
else
    error(['Please specify cfg.refdimension with ''braintime'' or ''clocktime''.'...
    ' See help bt_quantify or toolbox documentation for more details']);
end

% Perform analysis over over tgm, its autocorrelation map (ac), or its diagonal?
if (strcmp(maptype,'tgm') && strcmp(config.method,'tgm')) || strcmp(maptype,'diag')
    mp = mv_perf;
elseif (strcmp(maptype,'tgm') && strcmp(config.method,'ac')) 
    mp = autocorr2d(mv_perf);
end

[PS,f] = bt_fft(mp,powspecrange,timevec);
pspec_emp = PS;

figopt = bt_defaultval(config,'figure','yes');

% Plotting tgm diagonal (no time generalization)
if strcmp(figopt,'yes') && strcmp(maptype,'diag')
    % Plot Diagonal
    figure; bt_figure('halfwide');
    subplot(2,1,1)
    xvec = linspace(timevec(1),timevec(end),numel(mp));
    plot(xvec,mp,'LineWidth',3,'Color',bt_colorscheme('diagonal'));
    
    if isempty(bt_plotparams('diag_ylim'))~=1 % Change y-lim
    ylim(bt_plotparams('diag_ylim'));
    end
    
    try % For old Matlab versions
    yline(0.5,'LineWidth',1.5,'Color',[0.6 0.6 0.6]);
    catch
    vline(0.5); 
    vline(0.5); 
    end
    
    % xlabel is dependent on refdimension
    if strcmp(refdimension.dim,'braintime')
        xlabel('time (cycles)');
    elseif strcmp(refdimension.dim,'clocktime')
        xlabel('time (seconds)');
    end
    ylabel(['Performance (',mv_cfg.metric,')']);
    title('TGM diagonal (classifier time course)');
    
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    
    % Prepare second subplot
    hold on
    subplot(2,1,2);
    
    % Plotting tgms and ac maps
elseif strcmp(figopt,'yes') && strcmp(maptype,'tgm')
    % Plot tgm
    figure; bt_figure('halfwide');
    subplot(2,2,1);
    cfg_plot= [];
    cfg_plot.clim = bt_plotparams('TGM_clim');
    cfg_plot.x   = linspace(timevec(1),timevec(end),size(mp,1));
    cfg_plot.y   = linspace(timevec(1),timevec(end),size(mp,2));
    mv_plot_2D(cfg_plot, mv_perf);
    axis square
    cb = colorbar;
    colormap(bt_colorscheme('TGM'));
    title(cb,mv_cfg.metric);freezeColors;
    xlim([timevec(1) timevec(end)]);
    ylim([timevec(1) timevec(end)]);
    xticks(yticks) % make ticks the same on the two axes
    title(['Time Generalization Matrix'])
    if strcmp(refdimension.dim,'braintime')
        xlabel('Test data (cycles)')
        ylabel('Training data (cycles)')
    elseif strcmp(refdimension.dim,'clocktime')
        xlabel('Test data (seconds)')
        ylabel('Training data (seconds)')
    end
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    
    % Plot ac map
    mp2 = autocorr2d(mv_perf);
    
    % Has user specified a color limit?
    if isempty(bt_plotparams('AC_clim'))~=1
    clim = bt_plotparams('AC_clim');
    else 
    % Detect appropriate color range by zscoring
    mn_ac=median(mp2(:));
    sd_ac=std(mp2(:));
    ac_z=(mp2-mn_ac)/sd_ac;
    indx = (abs(ac_z)<5); %filter to include only z-scores under 5
    clim = [-max(mp2(indx)),+max(mp2(indx))]; %take the max number as clim for plotting
    end
    
    subplot(2,2,2);
    pcolor(timevec,timevec,mp2(1:numel(timevec),1:numel(timevec)));shading interp;title('Autocorrelation map');
    axis square
    colormap(bt_colorscheme('AC')); 
    caxis(clim);
    xticks(linspace(timevec(1),timevec(end),6)); % Create 11 steps
    yticks(xticks);
    xticklabels(linspace(-50,50,6)); % Indicate % shifted
    yticklabels(xticklabels);
    cb=colorbar;
    title(cb,'corr');freezeColors;
    xlabel('Map shifted by x%')
    ylabel('Map shifted by y%')
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    
    % Prepare second subplot
    hold on
    subplot(2,2,3:4)
end

% Now plot the power spectrum; which is done identically for all methods (only the subplot index is different)
p1 = plot(f,pspec_emp,'LineStyle','-','LineWidth',4,'Color',bt_colorscheme('per_ps_emp')); %Mean across 1st level perms
xlabel('Periodicity frequency (Hz)')
ylabel('Periodicity power')
xlim([f(1) f(end)]);

if strcmp(refdimension.dim,'braintime') %warp freq line is dependent on clock (warped freq) or brain time (1 hz)
    wfreq_i = nearest(f,1);
    p3 = line([f(wfreq_i) f(wfreq_i)], [0 max(pspec_emp)],'color',bt_colorscheme('warpingsignal'),'LineWidth',4); %Line at warped freq
    xlabel('Periodicity frequency (factor of warped freq)')
    p3.Color(4) = 0.85;
    
    % Add legend
    l1 = legend('Periodicity power','Warped frequency');
    l1.FontSize = bt_plotparams('FontSizeLegend');
end

if strcmp(maptype,'diag')
    title('Periodicity power spectrum of TGM diagonal (classifier time series)');
elseif strcmp(config.method,'ac')
    title('Periodicity power spectrum of TGM''s autocorrelation map');
elseif strcmp(config.method,'tgm')
    title('Periodicity power spectrum of TGM');
end

% Adapt font
set(gca,'FontName',bt_plotparams('FontName'));
set(gca,'FontSize',bt_plotparams('FontSize'));

%% Give saved warning messages
% Only given at the end, or else they are drown out by FieldTrip messages
if exist('warning_msg1','var') == 1
    warning(warning_msg1)
end
if exist('warning_msg2','var') == 1
    warning(warning_msg2)
end

%% Save basic info
quant.mv_results = mv_results;                    % MVPA Light results structure
quant.clabel = clabel;                            % Classification labels
quant.toi = toi;                                  % Start and end time of interest
quant.warpfreq = warpfreq;                        % Warped frequency (frequency of the warping signal)
quant.timevec = timevec;                          % Time vector (different for brain and clock time referencing)
quant.refdimension = refdimension;                % Reference dimension used
quant.periodicityfoi = powspecrange;              % Range of tested periodicity frequencies
quant.method = config.method;                     % Save whether analysis is done over TGM, AC map, or diag
quant.maptype = maptype;                          % MVPA output type (TGM or diagonal)
quant.pspec_emp = pspec_emp;                      % Periodicity power spectrum of empirical data
