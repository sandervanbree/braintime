function [fft_chans] = bt_analyzechannels(config,configFT,channels)
% Analyze the time frequency characteristics of channels (e.g.
% ICA components, virtual channels, LFP time series) in anticipation of
% brain time carrier extraction in bt_choosecarrier. Sorting is based 
% on the power in frequency bands of interest, and optionally the 
% correlation to a template topography (obtained with bt_templatetopo).

% Use:
% [fft_chans] = bt_analyzechannels(cfg,cfg.FT,channels)
%
% Input Arguments:
% config             % TOOLBOX configuration structure with cells:
%                    %
%   - time           % Start and end of the time window of interest
%   - warpfreqs      % Lowest and highest frequency of interest for the to
%                    % to be designated brain time oscillation (e.g.
%                    % 8 to 12 for alpha oscillations for attention)
%   - correct1f      % (Optional). 'yes' to apply 1/F correction to
%                    % visualized power spectrum (default: 'yes').
%   - ntopchans      % Number of best channels to be filtered from the full
%                    % amount.
%                    %
%   - cutmethod      % 'consistenttime': warp from start to end of window
%                    % of interest and no more. The upside of this method
%                    % is that the final brain time data is of the exact
%                    % intended duration. The downside is an artefact in
%                    % the first cycle caused by the dynamic time warping
%                    % algorithm.
%                    % 
%                    % 'cutartefact': warp half a second before and after
%                    % your time window of interest that is later cut. 
%                    % The upside is that the first cycle artefact is 
%                    % removed. The downside is variance across trials
%                    % in the brain time data's start and end time, as well
%                    % as the duration being slightly off to the window
%                    % of interest.
%                             
%   - sortmethod     % 'maxpow': sort channels according to their average
%                    % power in minfoi and maxfoi.
%                    %
%                    % 'templatetopo': loads a template topography created
%                    % using bt_templatetopo that is used to bias the
%                    % power sorting to the each channel's topography
%                    % match to the template topography.
%                    %
% configFT           % FIELDTRIP configuration structure with cells:
%                    %
%   - foi            % [min max]: Lowest and highest frequency of interest
%                    % analyzed by ft_freqanalysis. Choose frequencies
%                    % significantly below and above your warping frequency
%                    % of interest.
%   - method         % 'string': Method used by ft_freqanalysis ('mtmfft',
%                    % 'mtmconvol', 'wavelet', 'tfr', 'mvar'). See
%                    % ft_freqanalysis for details.
%   - pad            % (Optional). [number], 'nextpow2', or the default
%                    % 'maxperlen'. See ft_freqanalysis for details.
%   - padtype        % (Optional). 'string' of the padding type ('zero'
%                    % by default. See ft_freqanalysis for details.
%                    %
%                    % configFT also requires method-specific parameters
%                    % (e.g. 'tapsmofrq' for method 'mtmconvol'). Without
%                    % necessary parameters, FieldTrip will throw an error
%                    % and report on missing parameters. See
%                    % ft_freqanalysis for details.
%                    %
% channels           % FieldTrip channel data structure that contains
%                    % ICA components, virtual channels, or LFP time
%                    % series. One channel contains the to be designated
%                    % brain time carrier.
%                    %
% Output:            %
% fft_chans          % Data structure with: ranked channels, their 
%                    % time frequency information, and config details


%% Get information
sampledur = (channels.time{1}(2)-channels.time{1}(1));  % Duration of each sample
numchans = size(channels.trial{1},1);                   % Number of channels
minfoi = config.warpfreqs(1);                           % Lowest freq of interest
maxfoi = config.warpfreqs(2);                           % Highest freq of interest
mintime = config.time(1);                               % Minimum time of interest
maxtime = config.time(2);                               % Maximum time of interest
cfgFT = configFT;                                       % FieldTrip config structure
correct1f = config.correct1f;                           % Correct for 1/f in visualization of FFT option
ntopchans = config.ntopchans;                           % Number of top channels to be considered

% Amount of time dependent on cut method
if strcmp(config.cutmethod,'consistenttime')
    mintime_fft = mintime;
    maxtime_fft = maxtime;
elseif strcmp(config.cutmethod,'cutartefact')
    mintime_fft = mintime-0.5; 
    maxtime_fft = maxtime+0.5;
end

% Sanity check
if minfoi<cfgFT.foi(1) || maxfoi>cfgFT.foi(end)
error('Warping frequencies of interest detected outside specified FFT frequency range');    
end

%% Calculate FFT
cfgFT.toi       = mintime_fft:sampledur:maxtime_fft;
cfgFT.output    = 'fourier';
fspec           = ft_freqanalysis(cfgFT,channels);

fspec.powspctrm = abs(fspec.fourierspctrm);
fspec_old       = fspec;

if strcmp(correct1f,'yes')
fspec           = uni_subtract1f(fspec); % apply 1/F subtraction. This is just temporary to find the right channel
end

%% Find phase of chans in freq range of interest
% Find indices of interest
minfoi_ind = nearest(fspec.freq,minfoi); % index of lowest warping frequency
maxfoi_ind = nearest(fspec.freq,maxfoi); % index of highest warping frequency
foivec_ind = minfoi_ind:maxfoi_ind;      % freq range of interest vector
mintime_ind = nearest(fspec.time,mintime); % minimun time of interest index
maxtime_ind = nearest(fspec.time,maxtime); % maximun time of interest index

% Pre-allocate
powtf = zeros(size(fspec.powspctrm,3),maxtime_ind+1-mintime_ind,numchans);
pspec = zeros(size(fspec.powspctrm,3),numchans);
oscmaxfreq = zeros(numchans,1);
foi_ind = zeros(numchans,1);
phs = cell(numchans,1);

% For each channel, obtain time frequency information
for chan = 1:numchans
    powtf(:,:,chan)=squeeze(nanmean(fspec.powspctrm(:,chan,:,mintime_ind:maxtime_ind),1));
    pspec(:,chan)=squeeze(nanmean(powtf(:,:,chan),2)); % get power spectrum of channels
    [maxpow, maxpowind]= max(pspec(minfoi_ind:maxfoi_ind,chan)); % what's the highest power in the freq range of interest and its index in the freq range of interest vector?
    oscmaxfreq(chan) = (fspec.freq(foivec_ind(maxpowind))); % what's the highest power oscillation?
    [~, freqidx]=find(abs(oscmaxfreq(chan)-fspec.freq)==min(abs(oscmaxfreq(chan)-fspec.freq))); % what is the index of the highest power oscillation?
    foi_ind(chan) = freqidx;
    oscmax(chan,:) = [chan,oscmaxfreq(chan), foi_ind(chan),maxpow]; % matrix with channel number, freq, its index, and its power
    phs{chan}=squeeze(angle(fspec_old.fourierspctrm(:,chan,foi_ind(chan),:))); % fetch phase of highest power oscillation
end

%% Sort channels by one of two methods
% For maxpow this is sufficient:
chanrank = sortrows(oscmax,4,'descend'); % sort channels from highest to lowest power

% Template topography builds onto the variable chanrank
if strcmp(config.sortmethod,'templatetopo')
    try
        load temptopo
    catch
        error('Template topography not found. Use bt_templatetopo to create a template topography')
    end
    
    % create final template topography
    try
        topo_ind = contains(channels.topolabel,temptopo); % Find the indices of the chosen sensors in the channel's labels
    catch
        error('Channel labels in your template topography mismatch carrier labels')
    end
    topo = double(topo_ind)*5; %Put high value (5) only for chosen channels
    
    % Get the channel-template correlation
    corrord=chanrank;
    corrord(:,1)=1:numchans;
    
    % correlate template topography with channel topography
    for chan=1:numchans
        tc=corrcoef(topo,abs(channels.topo(:,chan)));
        corrord(chan,2)=tc(1,2);
    end
    corrord=sortrows(corrord,2,'descend'); % sort channels from highest to lowest correlation
    corrord=corrord(:,1);
    
    crank=numchans;
    for i=1:crank
        chanrank(i,5)=crank; % add descending values to the power order channels
        corrord(i,2)=crank; % add descending values to the correlation order channels
        crank=crank-1;
    end
    chanrank=sortrows(chanrank,1,'ascend'); % sort by channels number
    corrord=sortrows(corrord,1,'ascend');
    
    chanrank(:,6)=chanrank(:,5)+corrord(:,2); %get an index with the addition
    
    chanrank=sortrows(chanrank,6,'descend'); % sort channels by index
    chanrank =chanrank(:,1:4);
end
    
%% Take the 30 best channels (or specified number)
if isfield(config,'ntopchans')
    ntopchans = min(ntopchans,numchans);
else
    ntopchans = min(numchans,30);
end
chanrank = chanrank(1:ntopchans,:); % take the top numcarrier carriers

%% Save information
% Filter relevant frequency spectrum information
fspecinfo.freq = fspec.freq;
fspecinfo.time = fspec.time;

fft_chans{1} = chanrank;                                  % Time freq data of chosen channel
fft_chans{2} = [mintime_ind maxtime_ind];                 % Index of start and end time of interest (differs for cutartefact)
fft_chans{3} = [cfgFT.foi(1) cfgFT.foi(end)];             % Lowest and highest analyzed freq in FFT
fft_chans{4} = [minfoi maxfoi];                           % Lowest and highest possible warping frequency
fft_chans{5} = fspecinfo;                                 % FFT time and frequency vector
fft_chans{6} = powtf(:,:,chanrank(:,1));                  % Power spectrum of channels
fft_chans{7} = pspec(:,chanrank(:,1));                    % Power spectrum averaged across trials
fft_chans{8} = phs(chanrank(:,1));                        % Phase of all channels for all trials
fft_chans{9} = config.cutmethod;                          % Applied cutting method
end