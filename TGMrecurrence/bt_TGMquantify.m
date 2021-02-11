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
%   - mapmethod      % 'ac' (default): perform recurrence analysis over the
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
warpfreq = config.bt_warpeddata.freq;                 % Warped frequency (frequency of the carrier) 
duration = toi(2)-toi(1);                             % Duration of the time window of interest
mapmethod = config.mapmethod;                         % Perform the analysis over TGM or its autocorrelation map?

% Set up recurrence range over which stats will be applied
if isfield(config,'recurrencefoi')
    powspecrange = config.recurrencefoi(1):config.recurrencefoi(end);
else
    warning('No recurrence frequency range of interest entered under config.recurrencefoi; testing 1 to 30 Hz')
    powspecrange = 1:30;
end

if strcmp(config.refdimension,'braintime')
    powspecrange = powspecrange/warpfreq; % adjust tested range to warped frequency
    
    refdimension.value = duration*warpfreq; % normalize by cycles in the data
    refdimension.dim = 'braintime';
    timevec = config.bt_warpeddata.data.time{1};
elseif strcmp(config.refdimension,'clocktime')
    refdimension.value = duration; %normalize by seconds in the data
    refdimension.dim = 'clocktime';
    timevec = linspace(toi(1),toi(2),numel(config.bt_warpeddata.data.time{1}));
else
    error('Please specify cfg.refdimension with ''braintime'' or ''clocktime''. See help bt_TGMquantify or toolbox documentation for more details');
end

% Perform analysis over TGM or its autocorrelation map (AC)?
if strcmp(mapmethod,'tgm')
    mp = TGM;
elseif strcmp(mapmethod,'ac')
    mp = autocorr2d(TGM);
else
    warning('No config.mapmethod specified. Will perform analysis over TGM''s autocorrelation map (default)');
end   

% Run FFT over all rows and columns of the AC map
nrows = numel(mp(:,1));
ncols = numel(mp(:,2));

for row = 1:nrows % Perform FFT over rows
        
    if row == 1 % For the first row, perform a test analysis
        [~,f]=Powspek(mp(1,:),nrows/refdimension.value);
        l = nearest(f,powspecrange(1)); %minimum frequency to be tested
        h = nearest(f,powspecrange(end)); %maximum frequency to be tested
        ps_range = l:h; % this is the range of frequencies desired
    end
    
    % 1st dimension
    [PS,f]=Powspek(mp(row,:),nrows/refdimension.value);
    PS1(row,:) = PS(ps_range); % restrict do desired range
end

for col = 1:ncols % Perform FFT over columns
    
    if col == 1 % For the first column, perform a test analysis
        [~,f]=Powspek(mp(1,:),ncols/refdimension.value);
        l = nearest(f,powspecrange(1)); %minimum frequency to be tested
        h = nearest(f,powspecrange(end)); %maximum frequency to be tested
        ps_range = l:h; % this is the range of frequencies desired
    end
    
    % 2nd dimension
    [PS,f]=Powspek(mp(:,col),ncols/refdimension.value);
    PS2(col,:) = PS(ps_range); % restrict do desired range
end

f=f(ps_range); %filter frequency vector based on range of interest
avg_PS = mean(PS1,1)+mean(PS2,1); %Mean power spectra
pspec_emp = avg_PS;

if isfield(config,'figure')
    if strcmp(config.figure,'no')
        figopt = 0;
    else
        figopt = 1; %Default yes
    end
end

if figopt == 1
    % Plot TGM
    figure;
    subplot(2,2,1)
    cfg_plot= [];
    cfg_plot.x   = timevec;
    cfg_plot.y   = cfg_plot.x;
    mv_plot_2D(cfg_plot, TGM);
    cb = colorbar;
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
    % Detect appropriate color range by zscoring  
    mn_ac=median(mp(:));
    sd_ac=std(mp(:));
    ac_z=(mp-mn_ac)/sd_ac;
    indx = (abs(ac_z)<5); %filter to include only z-scores under 5
    clim = max(mp(indx)); %take the max number as clim for plotting
    
    subplot(2,2,2)
    pcolor(timevec,timevec,mp(1:numel(timevec),1:numel(timevec)));shading interp;title(['Autocorrelation map'])
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
    
    hold on
    
    subplot(2,2,3:4)
    p1 = plot(f,pspec_emp,'LineStyle','-','LineWidth',3,'Color','b'); %Mean across 1st level perms
    xlabel('Recurrence frequency (Hz)')
    ylabel('Mean power')
  
    if strcmp(refdimension.dim,'braintime') %warp freq line is dependent on clock (warped freq) or brain time (1 hz)
        p3 = line([1 1], [0 max(pspec_emp)],'color',[1 0 1],'LineWidth',4); %Line at warped freq
        xlabel('Recurrence frequency (factor of warped freq)')
    else
        p3 = line([warpfreq warpfreq], [0 max(pspec_emp)],'color',[1 0 1],'LineWidth',4); %Line at warped freq
        xlabel('Recurrence frequency')
    end
    p3.Color(4) = 0.45;
         
    if strcmp(mapmethod,'ac')
        title('Recurrence power spectrum of TGM''s autocorrelation map');
    else
        strcmp(mapmethod,'tgm')
        title('Recurrence power spectrum of TGM');
    end
    
    % Add legend
    legend('Recurrence power','Warped frequency');
    
    % Adapt font
    set(gca,'FontName','Arial')
    set(gca,'FontSize',16)
end

%% Save basic info
bt_TGMquant.toi = toi;                                  % Start and end time of interest
bt_TGMquant.warpfreq = warpfreq;                        % Warped frequency (frequency of the warping signal) 
bt_TGMquant.timevec = timevec;                          % Time vector (different for brain and clock time referencing)
bt_TGMquant.refdimension = refdimension;                % Reference dimension used
bt_TGMquant.recurrencefoi = powspecrange;               % Range of tested TGM recurrence frequencies
bt_TGMquant.TGM = TGM;                                  % Save the TGM itself
bt_TGMquant.mapmethod = mapmethod;                      % Save whether analysis is done over TGM or AC map
bt_TGMquant.pspec_emp = pspec_emp;                      % Recurrence power spectrum of empirical data
