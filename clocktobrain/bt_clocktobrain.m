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
%   - bt_srate       % Sampling rate of the brain time data.
%                    %
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
%                    % 'stationary': warp to a stationary sinusoid at the
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
phs = cell2mat(bt_source{2});                          % Phase of all frequencies in this warping source
warpsources = bt_source{3};                            % Warping source data
srcrank = bt_source{4};                                % Time freq data of selected warping signal
mintime_fft = bt_source{5}.time(1);                    % Start time of interest
maxtime_fft = bt_source{5}.time(end);                  % End time of interest
sr = bt_source{5}.time(2)-bt_source{5}.time(1);        % Sampling rate
cutmethod = bt_source{6};                              % Applied cutting method
warpfreq = srcrank(2);                                 % Warped frequency (frequency of the warping signal)
wavshape = bt_source{7};                               % Average waveshape of the data

warpmethod = bt_defaultval(config,'warpmethod','stationary');  % Set method for warping (default: stationary)

if strcmp(cutmethod,'cutartefact')                     % Depending on cutmethod, specify original time window of interest
    mintime = mintime_fft+0.5;
    maxtime = maxtime_fft-0.5;
else
    mintime = mintime_fft;
    maxtime = maxtime_fft;
end

mintime_ind = nearest(bt_source{5}.time,mintime);    % Index of start time of interest (differs for cutartefact)
maxtime_ind = nearest(bt_source{5}.time,maxtime);    % Index of end time of interest

% Set up sampling rate
if isfield(config,'bt_srate')
    phs_sr = config.bt_srate;
else
    phs_sr = 512; %Default sampling rate
end

%% Remove the component from original data if desired (default = yes)
cfg           = [];
cfg.component = src_oi;
if isfield(config,'removecomp')
    if strcmp(config.removecomp,'yes')
        data = ft_rejectcomponent (cfg, warpsources, data);
    end
else
    % if the removal option is not specified, remove by default
    if strcmp(warpsources.cfg.method,'runica')
        data = ft_rejectcomponent (cfg, warpsources, data);
    end
end

%% Cut out the time window of fft (from which the phase was extracted)
cfg        = [];
if strcmp(cutmethod,'cutartefact')
    cyclesample = round((1/warpfreq)*1/sr); % Calculate how many samples one cycle consists of
    cfg.toilim = [mintime_fft+0.5-(1/warpfreq) maxtime_fft-0.5+(1/warpfreq)]; % Cut to the time window of interest, plus one cycle
    phs = phs(:,mintime_ind-cyclesample:maxtime_ind+cyclesample); % Cut to the time window of interest, plus one cycle
elseif strcmp(cutmethod,'consistenttime')
    cfg.toilim = [mintime_fft maxtime_fft];
end
data       = ft_redefinetrial(cfg, data);

% Check whether phase and data are of the same length
if abs(length(data.time{1})-size(phs,2))>1
    warning('phase vector and data differ in length by more than 1 sample. This may mean the wrong parameters were used during bt_analyzesources');
end
if size(phs,2) > length(data.time{1})
    phs=phs(:,1:length(data.time{1}));
elseif length(data.time{1}) > size(phs,2)
    cfg       = [];
    cfg.toilim = [data.time{1}(1,1) data.time{1}(1,size(phs,2))];
    data      = ft_redefinetrial(cfg, data);
end

%% Warp the phase of the warping signal to the phase of a stationary oscillation
bt_data=data;
nsec=bt_data.time{1}(end)-bt_data.time{1}(1);                 % number of seconds in the data
Ncycles_pre=warpfreq*nsec;                                    % number of cycles * seconds
cycledur=round(phs_sr*nsec/Ncycles_pre);                      % samples for cycle
tempsr=Ncycles_pre*cycledur/nsec;

if strcmp(warpmethod,'stationary')                            % warp using stationary sinusoid
    tempphs=linspace(-pi,(2*pi*Ncycles_pre)-pi,tempsr*nsec);  % set up phase bins for unwrapped phase (angular frequency)
    
elseif strcmp(warpmethod,'waveshape')                         % warp using average waveshape
    waveshape = smoothdata(wavshape,'gaussian',100);          % smooth substantially
    [~,trgh] = findpeaks(-waveshape);                         % find troughs
    
    if numel(trgh)~=2                                         % if there are not 2 troughs, this likely
        error(['The waveshape of the warping signal is',...   % means the data is too noisy
            ' too noisy. Please select cfg.method =',...
            ' ''stationary'' and try again.']);
    end
    
    waveshape_cut = waveshape(trgh(1)+1:trgh(2));             % cut waveshape to one cycle (-cos)
    tempsig = repmat(waveshape_cut,[1 floor(Ncycles_pre)]);   % repeat wave Ncycles_pre times (round down)
    misscycs = Ncycles_pre-floor(Ncycles_pre);                % check if rounding down lost us anything
    extrasamps = round(numel(waveshape_cut)*misscycs);        % how many samples of the cut cycle were lost
    
    if misscycs~=0                                            % append those to the template signal
    tempsig = [tempsig, waveshape_cut(1:extrasamps)];
    end
    
    tempphs = angle(hilbert(tempsig));                        % get the phase using Hilbert transform
    tempphs = unwrap(tempphs);                                % unwrap the phase
    tempphs = imresize(tempphs,[1 tempsr*nsec]);              % resize to the desired length
end

timephs=linspace(0,Ncycles_pre,phs_sr*nsec);                  % time vector of the unwrapper phase

for nt=1:size(phs,1)
    tmpphstrl=unwrap((phs(nt,:)));
    % Warp phase of single trial onto template phase
    [~,ix,iy] = dtw(tmpphstrl,tempphs);
    
    % How long is each cycle?
    cycles=zeros(1,length(iy));
    cycles(1,end)=500;
    for tp=1:length(iy)-1
        if rem(iy(tp),cycledur)==0
            if iy(tp)~=iy(tp+1)
                cycles(tp)=500; % Label the start of each cycle with an arbitrary '500'
            end
        end
    end
    [~, c]=find(cycles==500); % Ffind the indices of the starts
    
    % Create an equal number of samples per cycle
    % Code for the first cycle
    cyl=bt_data.trial{1,nt}(:,ix(1:c(1)));
    tmpcy=imresize(cyl,[size(bt_data.label,1) cycledur]);
    tmptrl(:,1:cycledur)=tmpcy;
    
    % A loop for the remaining cycles
    for cy=2:Ncycles_pre
        cyl=bt_data.trial{1,nt}(:,ix(c(cy-1)+1:c(cy)));
        tmpcy=imresize(cyl,[size(bt_data.label,1) cycledur]);
        tmptrl(:,(cy-1)*cycledur+1:cy*cycledur)=tmpcy;
    end
    
    % Create brain time warped trials by resizing to original length
    bt_data.trial{1,nt}=imresize(tmptrl,[size(tmptrl,1) numel(timephs)]);
    bt_data.time{1,nt}=timephs;
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

%% Reformat data structure and include basic info
bt_warpeddata.data = bt_data;                                     % Brain time warped data
bt_warpeddata.toi = [min_t max_t];                                % Start and end time of interest
bt_warpeddata.freq = warpfreq;                                    % Warped frequency (frequency of the warping signal)
bt_warpeddata.clabel = bt_data.trialinfo;                         % Classification labels
bt_warpeddata.method = warpmethod;                                    % Warping method