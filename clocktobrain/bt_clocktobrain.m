function [bt_warpeddata] = bt_clocktobrain(config, data, bt_source)
% Warp clock to brain time. The phase of the warping signal in the selected
% warping source is dynamically time warped to the phase of a stationary
% oscillation. This warping path is used to resample the clock time data,
% stretching at specific moments where the warping signal falls out of tune
% with the stationary oscillation. Finally, the data is resized to its
% original length. This procedure attempts to account for variation in each
% trial's starting phase as well as frequency drift within trials.
%
% Use:
% [bt_struc] = bt_clocktobrain(cfg,data,bt_source)
%
% Input Arguments:
% config
%   - removecomp     % 'yes' (default): if the warping source is an ICA
%                    % component, remove it from the brain time data. When
%                    % analyzing brain time data using your own analysis
%                    % pipeline you may wish to remove the component to
%                    % avoid circularity. See the brain time toolbox paper
%                    % for more details on issues of circularity.
%                    %
%                    % 'no': keeps component in the brain time data.
%                    %
%   - method         % warping method: the desired template oscillation
%                    % for warping.
%                    %
%                    % 'sinusoid': warp to a statonary sinusoid at the
%                    % warping frequency (default).
%                    %
%                    % 'waveshape': warp to the average waveshape of the
%                    % warping signal. Recommended for asymmetric warping
%                    % signals.
%                    %
% data               % Clock time data structure consisting of both
%                    % classes.
%                    %
% bt_source          % Data structure obtained from bt_selectsource.
%                    % Includes: Warping signal's time frequency
%                    % information, and config details saved for later
%                    % retrieval.
%                    %
% Output:            %
% bt_warpeddata      % Data structure with: the warping signal's frequency,
%                    % brain time data, its time frequency information, and
%                    % config details saved for later retrieval.

%% Get basic info
src_oi = bt_source{1};                                 % Warping source which contains the warping signal
FFT_phs = cell2mat(bt_source{2});                      % Phase of all frequencies in this warping source
GED_phs = bt_source{9};                                % Phase of warping signal as estimated using GED
warpsources = bt_source{3};                            % Warping source data
srcrank = bt_source{4};                                % Time freq data of selected warping signal
mintime_fft = bt_source{5}.time(1);                    % Start time of interest
maxtime_fft = bt_source{5}.time(end);                  % End time of interest
sr = bt_source{5}.time(2)-bt_source{5}.time(1);        % Sampling rate
cutmethod = bt_source{6};                              % Applied cutting method
warpfreq = srcrank(2);                                 % Warped frequency (frequency of the warping signal)
wvshape = bt_source{7};                                % Average waveshape of the data

warpmethod = bt_defaultval(config,'warpmethod','sinusoid');    % Set method for warping (default: sinusoid)
phasemethod = bt_defaultval(config,'phasemethod','FFT');       % Set phase estimation method used for warping (default: FFT)
visualcheck = bt_defaultval(config,'visualcheck','off');       % Show warping path and phases for three example trials

if strcmp(cutmethod,'cutartefact')                     % Depending on cutmethod, specify original time window of interest
    mintime = mintime_fft+0.5;
    maxtime = maxtime_fft-0.5;
else
    mintime = mintime_fft;
    maxtime = maxtime_fft;
end

mintime_ind = nearest(bt_source{5}.time,mintime);    % Index of start time of interest (differs for cutartefact)
maxtime_ind = nearest(bt_source{5}.time,maxtime);    % Index of end time of interest

% Transform data structure to raw
if isstruct(data.trial) == 0
data = ft_checkdata(data,'datatype','raw');
warning('Timelocked data structure detected, so it was converted to raw');
end

% Set sampling rate of brain time to clock time sampling rate
time   = data.time{1};
phs_sr = round(1/(time(2)-time(1)));

%% Remove the component from original data if desired (default = yes)
circwarning   = 1;        % Warn the user for potential circularity by default

cfg           = [];
cfg.component = src_oi;
if isfield(config,'removecomp')
    if strcmp(config.removecomp,'yes')
        if isfield(warpsources,'cfg') && isfield(warpsources.cfg,'method')
            if strcmp(warpsources.cfg.method,'runica')
                data = ft_rejectcomponent (cfg, warpsources, data);
                circwarning = 0; % Component removed, so no need for a warning
            else
                error(['It is not possible to remove the warping component, because'...
                    ' the warping sources are not ICA components.']);
            end
        else
            error(['It is not possible to remove the warping component, because'...
                ' the warping sources are not ICA components.']);
        end
    end
else
    % if the removal option is not specified, remove by default
    if isfield(warpsources,'cfg') && isfield(warpsources.cfg,'method') && strcmp(warpsources.cfg.method,'runica')
        data = ft_rejectcomponent (cfg, warpsources, data);
        circwarning = 0;
    end
end

%% Cut out the time window of phase estimation used during bt_analyzesources
cfg        = [];
if strcmp(cutmethod,'cutartefact')
    cyclesample = round((1/warpfreq)*1/sr); % Calculate how many samples one cycle consists of
    cfg.toilim = [mintime_fft+0.5-(1/warpfreq) maxtime_fft-0.5+(1/warpfreq)]; % Cut to the time window of interest, plus one cycle
    FFT_phs = FFT_phs(:,mintime_ind-cyclesample:maxtime_ind+cyclesample); % Cut to the time window of interest, plus one cycle
    GED_phs = GED_phs(:,mintime_ind-cyclesample:maxtime_ind+cyclesample); % Do the same for GED phase
elseif strcmp(cutmethod,'consistenttime')
    cfg.toilim = [mintime_fft maxtime_fft];
end
data       = ft_redefinetrial(cfg, data);


%% Opt for FFT or GED phase and adapt data
if strcmp(phasemethod,'FFT')
    phs = FFT_phs;
elseif strcmp(phasemethod,'GED')
    
    % Sanity check whether the two phase vectors are within 5% of each other's length
    if abs(size(FFT_phs,2)-size(GED_phs,2)) > 0.05*max(size(FFT_phs,2),size(GED_phs,2))
        error(['The length of the GED estimated phase substantially differs from '...
            'the length of the FFT estimated phase. Please change to config.phasemethod '...
            '= ''FFT'' or change the time window tested during bt_analyzechannels.'])
    else
        phs = GED_phs;
    end
end

% Check whether phase and data are of the same length
if abs(length(time)-size(phs,2))>1
    warning(['The phase and data length in length by more than 1 samples. This'...
        ' may mean the wrong parameters were used during bt_analyzesources.']);
elseif abs(length(time)-size(phs,2))>10
    error(['The phase and data length in length by more than 10 samples. This'...
        ' may mean the wrong parameters were used during bt_analyzesources.']);
end

% Do some additional slicing to fix slight differences
if size(phs,2) > length(time)
    phs=phs(:,1:length(time));
elseif length(time) > size(phs,2)
    cfg       = [];
    cfg.toilim = [time(1) time(end)];
    data      = ft_redefinetrial(cfg, data);
end

%% Warp the phase of the warping signal to the phase of a stationary oscillation
bt_data     = data;
nsec        = bt_data.time{1}(end)-bt_data.time{1}(1);        % number of seconds in the data
Ncycles     = warpfreq*nsec;                                  % number of cycles
cycledur    = round(phs_sr*nsec/Ncycles);                     % samples for cycle
tempsr      = Ncycles*cycledur/nsec;
timephs     = linspace(0,Ncycles,phs_sr*nsec);                % time vector of the unwrapper phase

if strcmp(warpmethod,'sinusoid')                              % warp using stationary sinusoid
    ct_phs = linspace(-pi,(2*pi*Ncycles)-pi,tempsr*nsec);     % set up phase bins for unwrapped phase (angular frequency)
    
elseif strcmp(warpmethod,'waveshape')                         % warp using average waveshape
    wvshape_sm = smoothdata(wvshape,'gaussian',25);           % smooth substantially - 25 bins seems the sweetspot for 201 bins
    [~,trgh] = findpeaks(-wvshape_sm);                        % find troughs
    
    if numel(trgh)~=2                                         % if there are not 2 troughs, this likely
        error(['The waveshape of the warping signal is',...   % means the data is too noisy
            ' too noisy. Please select cfg.method =',...
            ' ''sinusoid'' and try again.']);
    end
    
    waveshape_cut = wvshape_sm(trgh(1)+1:trgh(2));            % cut waveshape to one cycle (-cos)
    tempsig = repmat(waveshape_cut,[1 floor(Ncycles)]);       % repeat wave Ncycles times (round down)
    misscycs = Ncycles-floor(Ncycles);                        % check if rounding down lost us anything
    extrasamps = round(numel(waveshape_cut)*misscycs);        % how many samples of the cut cycle were lost
    
    if misscycs~=0                                            % append those to the template signal
        tempsig = [tempsig, waveshape_cut(1:extrasamps)];
    end
    
    tempsig_hb = angle(hilbert(tempsig));                     % get the phase using Hilbert transform
    ct_phs = unwrap(tempsig_hb);                              % unwrap the phase
    ct_phs = imresize(ct_phs,[1 tempsr*nsec]);                % resize to the desired length
    
    if strcmp(visualcheck,'on') || strcmp(visualcheck,'yes')  % perform visual check
        figure; hold on; bt_figure('halflong');
        
        subplot(3,1,1);
        tempsig_check = imresize(tempsig,[1 numel(timephs)]);                % Resize to timephs for plotting
        plot(timephs,tempsig_check,'LineWidth',3,'Color',[0.8 0.1 0.1]);     % Smoothed repeated waveshape
        title('Smoothed repeated waveshape');
        ylabel('Amplitude');
        
        subplot(3,1,2);
        tempsig_hb_check = imresize(tempsig_hb,[1 numel(timephs)]);
        plot(timephs,tempsig_hb_check,'LineWidth',3,'Color',[0.8 0.1 0.1]);  % Phase (Hilbert transform)
        title('Phase (Hilbert transform)');
        ylabel('Phase (-2*pi to 2*pi)');
        
        subplot(3,1,3);
        tempphs_check = imresize(ct_phs,[1 numel(timephs)]); 
        plot(timephs,tempphs_check,'LineWidth',3,'Color',[0.8 0.1 0.1]);     % Unwrapped phase
        title('Unwrapped phase');
        xlabel('Time (cycles or seconds)');
        ylabel('Unwrapped phase');
        
        % Adapt font
        set(findobj(gcf,'type','axes'),'FontName',bt_plotparams('FontName'),'FontSize',bt_plotparams('FontSize'));
    end
end

for nt=1:size(phs,1)
    bt_phs = unwrap((phs(nt,:)));
    % Warp phase of single trial onto template phase
    [~,ix,iy] = dtw(bt_phs,ct_phs);
    
    % Find the start of each cycle in the clock to brain time path
    cycles        = false(1,length(iy));
    cycles(1,end) = true;
    for p = 1:length(iy)-1            % Loop through path
        if rem(iy(p),cycledur) == 0   % If we have reached cycledur
            if iy(p) ~= iy(p+1)       % And the next bin is not a repetition
                cycles(p) = true;     % Label the start of the cycle
            end
        end
    end
   [~, c]=find(cycles == true);      % Find the indices of cycle starts
    
    % Resample the first cycle's data based on the warping path from brain
    % to clock time
    % First cycle:
    cyl   = bt_data.trial{1,nt}(:,ix(1:c(1)));
    tmpcy = imresize(cyl,[size(bt_data.label,1) cycledur]);  % Resize to cycle duration
    tmptrl(:,1:cycledur) = tmpcy;                            % Replace data
    
    % Remaining cycles:
    for cy    = 2:round(Ncycles)
        cyl   = bt_data.trial{1,nt}(:,ix(c(cy-1)+1:c(cy)));
        tmpcy = imresize(cyl,[size(bt_data.label,1) cycledur]);
        tmptrl(:,(cy-1)*cycledur+1:cy*cycledur) = tmpcy;
    end
    
    warpedtrial = imresize(tmptrl,[size(tmptrl,1) numel(timephs)]);
    
    % Perform visual check on warping path, phase, and warped data (first two trials)
    if (strcmp(visualcheck,'on') || strcmp(visualcheck,'yes')) && nt <= 2
        figure; hold on; bt_figure(0);
        dtw(bt_phs,ct_phs);     % warping path plot
        
        title(['Alignment after warping (trial ',num2str(nt), ')']);
        legend('Brain time phase',['Clock time phase (method: ',warpmethod,')'],'Location','northwest');
        xlabel('Data sample');
        ylabel('Unwrapped phase');
        
        % Adapt font
        set(findobj(gcf,'type','axes'),'FontName',bt_plotparams('FontName'),'FontSize',bt_plotparams('FontSize'));
        set(findobj(gcf,'type','line'),'LineWidth',1.5);

        figure; hold on; bt_figure(0);
        plot_prewarp = squeeze(mean(bt_data.trial{1,nt},1));  % Data before warping
        plot_postwarp = squeeze(mean(warpedtrial,1));         % Data after warping
        
        if size(plot_prewarp,2) > size(plot_postwarp,2)       % Crop to same length
            plot_prewarp = plot_prewarp(1:size(plot_postwarp,1),1:size(plot_postwarp,2));
        elseif size(plot_prewarp,2) < size(plot_postwarp,2)
            plot_postwarp = plot_postwarp(1:size(plot_prewarp,1),1:size(plot_prewarp,2));
        end
        
        plot(plot_prewarp,'LineWidth',3,'Color',[0.8 0.1 0.8]);
        plot(plot_postwarp,'LineWidth',3,'Color',[0.1 0.8 0.8]);
        
        if strcmp(cutmethod,'cutartefact')
        txt1 = '\leftarrow data cut from here';
        txt2 = 'and cut to here \rightarrow';
        text(cycledur,max(plot_postwarp),txt1,'FontSize',14,'HorizontalAlignment','left','FontWeight','bold');
        text(length(plot_postwarp)-cycledur,max(plot_postwarp),txt2,'FontSize',14,'HorizontalAlignment','right','FontWeight','bold');
        end
        
        title(['Average across channels (trial ',num2str(nt), ')']);
        legend('Pre-warp data','Post-warp data','Location','best');
        xlabel('Data sample');
        ylabel('Amplitude');
        
    end
    
    
    % Create brain time warped trials by resizing to original length
    bt_data.trial{1,nt} = warpedtrial;
    bt_data.time{1,nt}  = timephs;
end

% If cutmethod is cutartefact, slice out the right time window and adjust time vector
if strcmp(cutmethod,'cutartefact')
    % Correct time window of interest to true time
    min_t = mintime_fft+0.5;
    max_t = maxtime_fft-0.5;
    
    startind = nearest(bt_data.time{1},1); % Find index of first cycle in window of interest
    endind = nearest(bt_data.time{1},warpfreq*(max_t-min_t)+1); %Find index of last cycle in window of interest
    
    cfg         = [];
    cfg.latency = [bt_data.time{1}(startind) bt_data.time{1}(endind)]; % Cut to time window of interest
    bt_data  = ft_selectdata(cfg,bt_data);
    bt_data.trialinfo = data.trialinfo;
    
    % Correct the cycles time vector
    for trl = 1:numel(bt_data.trial)
        bt_data.time{trl} = bt_data.time{trl}-1; % Cycles vector is off by 1
    end
else
    min_t = mintime_fft;
    max_t = maxtime_fft;
end

% Warn user for circularity
if circwarning == 1
    if isfield(warpsources,'cfg') && isfield(warpsources.cfg,'method') && strcmp(warpsources.cfg.method,'runica')
        warning(['The warping sources look like ICA components obtained from the'...
            ' clock time data. When using brain time transformed data for analyses outside of'...
            ' the toolbox, it is strongly recommended to remove the component'...
            ' by specifying cfg.removecomp = ''yes''.'])
    else
        warning(['When using brain time transformed data for analyses outside of the toolbox,'...
            ' please ensure that the clock time data and warping sources used in'...
            ' the toolbox were sufficiently independent. This is important to'...
            ' avoid circularity issues.'])
    end
end

%% Reformat data structure and include basic info
bt_warpeddata.data = bt_data;                                     % Brain time warped data
bt_warpeddata.toi = [min_t max_t];                                % Start and end time of interest
bt_warpeddata.freq = warpfreq;                                    % Warped frequency (frequency of the warping signal)
bt_warpeddata.clabel = bt_data.trialinfo;                         % Classification labels
bt_warpeddata.warpmethod = warpmethod;                            % Warping method (waveshape or stationary sinusoid)
bt_warpeddata.phasemethod = phasemethod;                          % Phase estimation method (FFT or GED)
