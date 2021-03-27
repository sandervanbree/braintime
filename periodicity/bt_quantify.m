function [bt_TGMquant] = bt_TGMquantify(config, TGM)
% Quantify the degree of recurrence in the time generalization matrix
% (TGM). "Recurrence" refers to fluctuating patterns of classification
% performance, which reflect fluctuations of the studied cognitive process.
% bt_TGMquantify displays the TGM, its autocorrelation map, and peforms
% a fast Fourier transform (FFT) over all rows and columns of the TGM.
% This quantifies recurrence in the TGM, and is necessary for first level
% statistics (bt_TGMstatslevel1).
%
% Use:
% [bt_TGMquant] = bt_TGMquantify(cfg, TGM)
%
% Input Arguments:
% config
%   - bt_warpeddata  % brain time data structure as obtained by
%                    % bt_clocktobrain.
%                    %
%   - MVPAcfg        % Load the cfg file used to generate the TGM. This
%                    % provides information to the toolbox about which
%                    % classification parameters were used, and it ensures
%                    % later functions create permuted TGMs using the same
%                    % classification parameters.
%                    %
%   - refdimension   % 'clocktime': find TGM recurrence as a function of
%                    % clock time seconds in the brain time data.
%                    %
%                    % 'braintime': find TGM recurrence as a function of
%                    % cycles of the warped frequency in the brain time
%                    & data.
%                    %
%   - recurrencefoi  % '[min max]': Range of TGM recurrence frequencies to
%                    % be quantified and statistically tested.
%                    %
%                    %
%   - mapmethod      % 'diag' (default): perform recurrence analysis over
%                    % the diagonal of the TGM. This ignores cross-temporal
%                    % recurrence.
%                    %
%                    % 'ac': perform recurrence analysis over the
%                    % autocorrelation map of the TGM. This accentuates
%                    % the primary recurrence frequency in the TGM, but may
%                    % drown out other frequencies.
%                    %
%                    % 'tgm': perform recurrence analysis over the TGM
%                    % itself. This method is more sensitive to TGMs with
%                    % multiple recurrence frequencies.
%                    %
%   - figure         % 'yes' (default) or 'no': display TGM and its AC map.
%                    %
% TGM                % TGM obtained by MVPA Light's mv_classify_timextime
%                    %
% Output:            %
% bt_quantTGM        % Data structure with: TGM or AC map, recurrence
%                    % power spectrum, and config details saved for later
%                    % retrieval.

%% Get basic info
toi = config.bt_warpeddata.toi;                       % Start and end time of interest
warpfreq = config.bt_warpeddata.freq;                 % Warped requency (frequency of the carrier)
duration = toi(2)-toi(1);                             % Duration of the time window of interest
mapmethod = config.mapmethod;                         % Perform the analysis over TGM or its autocorrelation map?
MVPAcfg = config.MVPAcfg;                             % MVPA Light configuration structured used to obtain TGM
refdimension.dim = config.refdimension;               % Extract reference dimension (clock or brain time)

% Set up recurrence range over which stats will be applied
if isfield(config,'recurrencefoi')
    powspecrange = config.recurrencefoi(1):config.recurrencefoi(end);
else
    warning('No recurrence frequency range of interest entered under config.recurrencefoi; testing 1 to 30 Hz')
    powspecrange = 1:30;
end

if strcmp(refdimension.dim,'braintime') % Test if range is OK
    if powspecrange(1)/warpfreq > 0.5
        answer = questdlg('The lowest frequency of the specified recurrence range is too high to test recurrence at half the warped frequency, which is an important harmonic. Please consider increasing it before continuing.', ...
            'Warning', ...
            'Continue','Cancel','Continue');
        
        if strcmp(answer,'Cancel')
            error('Please change cfg.recurrencefoi')
        end
    elseif powspecrange(end)/warpfreq < 2
        
        answer = questdlg('The lowest frequency of the specified recurrence range is too low to test recurrence at double the warped frequency, which is an important harmonic. Please consider increasing it before continuing.', ...
            'Warning', ...
            'Continue','Cancel','Continue');
        
        if strcmp(answer,'Cancel')
            error('Please change cfg.recurrencefoi')
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
    timevec = linspace(toi(1),toi(2),size(TGM,1));
else
    error('Please specify cfg.refdimension with ''braintime'' or ''clocktime''. See help bt_TGMquantify or toolbox documentation for more details');
end

% Perform analysis over TGM, its autocorrelation map (AC), or its diagonal?
if strcmp(mapmethod,'tgm')
    mp = TGM;
elseif strcmp(mapmethod,'ac')
    mp = autocorr2d(TGM);
elseif strcmp(mapmethod,'diag')
    mp = diag(TGM)';
else
    warning('No valid config.mapmethod specified. Will perform analysis over TGM''s diagonal; i.e. the classifier timecourse without time generalization (default)');
    mapmethod = 'diag';
    mp = diag(TGM)';
end

[PS,f] = fftTGM(mp,powspecrange,timevec);
pspec_emp = PS;

if isfield(config,'figure')
    if strcmp(config.figure,'no')
        figopt = 0;
    else
        figopt = 1; %Default yes
    end
end

% Plotting TGM diagonal (no time generalization)
if figopt == 1 && strcmp(mapmethod,'diag')
    % Plot Diagonal
    figure;
    subplot(2,1,1)
    xvec = linspace(timevec(1),timevec(end),numel(mp));
    plot(xvec,mp,'LineWidth',3,'Color',[0 0 0]);
    
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
    ylabel(['Performance (',MVPAcfg.metric,')']);
    title('TGM Diagonal (classifier time course)');
    
    % Adapt font
    set(gca,'FontName','Arial');
    set(gca,'FontSize',16);
    
    % Prepare second subplot
    hold on
    subplot(2,1,2);
    
    % Plotting TGMs and AC maps
elseif figopt == 1 && strcmp(mapmethod,'tgm') || strcmp(mapmethod,'ac')
    % Plot TGM
    figure;
    subplot(2,2,1)
    cfg_plot= [];
    cfg_plot.x   = linspace(timevec(1),timevec(end),size(mp,1));
    cfg_plot.y   = linspace(timevec(1),timevec(end),size(mp,2));
    mv_plot_2D(cfg_plot, TGM);
    cb = colorbar;
    colormap(flipud(brewermap([],'RdBu')));
    title(cb,'perf')
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
    set(gca,'FontName','Arial')
    set(gca,'FontSize',16)
    
    % Plot AC map
    mp2 = autocorr2d(TGM);
    
    % Detect appropriate color range by zscoring
    mn_ac=median(mp2(:));
    sd_ac=std(mp2(:));
    ac_z=(mp2-mn_ac)/sd_ac;
    indx = (abs(ac_z)<5); %filter to include only z-scores under 5
    clim = max(mp2(indx)); %take the max number as clim for plotting
    
    subplot(2,2,2)
    pcolor(timevec,timevec,mp2(1:numel(timevec),1:numel(timevec)));shading interp;title('Autocorrelation map');
    colormap(flipud(brewermap([],'RdBu')));
    caxis([-clim +clim])
    xticks(linspace(timevec(1),timevec(end),6)); % Create 11 steps
    yticks(xticks);
    xticklabels(linspace(-50,50,6)); % Indicate % shifted
    yticklabels(xticklabels);
    cb=colorbar;
    title(cb,'corr')
    xlabel('Map shifted by x%')
    ylabel('Map shifted by y%')
    % Adapt font
    set(gca,'FontName','Arial')
    set(gca,'FontSize',16)
    
    % Prepare second subplot
    hold on
    subplot(2,2,3:4)
end

% Now plot the power spectrum; which is done identically for all mapmethods (only the subplot index is different)
p1 = plot(f,pspec_emp,'LineStyle','-','LineWidth',3,'Color','b'); %Mean across 1st level perms
xlabel('Recurrence frequency (Hz)')
ylabel('Mean power')

if strcmp(refdimension.dim,'braintime') %warp freq line is dependent on clock (warped freq) or brain time (1 hz)
    p3 = line([1 1], [0 max(pspec_emp)],'color',[1 0 1],'LineWidth',4); %Line at warped freq
    xlabel('Recurrence frequency (factor of warped freq)')
    p3.Color(4) = 0.45;
    
    % Add legend
    legend('Recurrence power','Warped frequency');
end

if strcmp(mapmethod,'diag')
    title('Recurrence power spectrum of TGM diagonal (classifier time series)');
elseif strcmp(mapmethod,'ac')
    title('Recurrence power spectrum of TGM''s autocorrelation map');
elseif strcmp(mapmethod,'tgm')
    title('Recurrence power spectrum of TGM');
end

% Adapt font
set(gca,'FontName','Arial')
set(gca,'FontSize',16)

%% Save basic info
bt_TGMquant.MVPAcfg = MVPAcfg;                          % MVPA Light configuration structured used to obtain TGM
bt_TGMquant.toi = toi;                                  % Start and end time of interest
bt_TGMquant.warpfreq = warpfreq;                        % Warped frequency (frequency of the warping signal)
bt_TGMquant.timevec = timevec;                          % Time vector (different for brain and clock time referencing)
bt_TGMquant.refdimension = refdimension;                % Reference dimension used
bt_TGMquant.recurrencefoi = powspecrange;               % Range of tested recurrence frequencies
bt_TGMquant.TGM = TGM;                                  % Save the TGM itself
bt_TGMquant.mapmethod = mapmethod;                      % Save whether analysis is done over TGM, AC map, or diag
bt_TGMquant.pspec_emp = pspec_emp;                      % Recurrence power spectrum of empirical data
