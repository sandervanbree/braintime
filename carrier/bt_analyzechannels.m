function [fft_chans] = bt_analyzechannels(config, channels)
% Analyze the time frequency characteristics of channels (e.g.
% ICA components, virtual channels, LFP time series) in anticipation of
% brain time carrier extraction in bt_choosecarrier. Sorting is based 
% on the power in frequency bands of interest, and optionally the 
% correlation to a template topography (obtained with bt_templatetopo).

% Use:
% [fft_chans] = bt_analyzechannels(cfg, channels)
%
% Input Arguments:
% config
%   - time           % Start and end of the time window of interest
%   - fft            % Lowest and highest frequency analyzed in the
%                    % channels 
%   - foi            % Lowest and highest frequency of interest for the to
%                    % to be designated brain time oscillation (e.g.
%                    % 8 to 12 for alpha oscillations for attention)
%   - Ntop           % Number of best channels to be filtered from the full
%                    % amount.
%                    
%   - cutmethod      % 'consistenttime': warp from mintime to maxtime.
%                    % The upside of this method is that the final
%                    % brain time data is equally long and of the same
%                    % duration as intended (mintime to maxtime). The
%                    % downside is an artefact in the first cycle
%                    % caused by the warping requiring a repetition of data
%                    % bins for alignment to the template oscillation.
%                    % 
%                    % 'cutartefact': warp from mintime-one cycle to
%                    % maxtime+one cycle of the minimum frequency
%                    % of interest, later cut to mintime maxtime.
%                    % The upside of this method is that the final brain
%                    % time data has no artefact in the first cycle. The
%                    % downside is that there is variance across trials
%                    % in the brain time's start and end time, as well as
%                    % the duration.
%                             
%   - sortmethod     % 'maxpow': sort channels according to their average
%                    % power in minfoi and maxfoi.
%                    %
%                    % 'templatetopo': loads a template topography created
%                    % using bt_templatetopo that is used to bias the
%                    % power sorting to the each channel's topography
%                    % match to the template topography.
%                    %
% channels           % FieldTrip channel data structure that contains
%                    % ICA components, virtual channels, or LFP time
%                    % series. One channel contains the to be designated
%                    % brain time carrier.
%                    %
% Output:            %
% fft_chans          % Data structure with: ranked channels, their 
%                    % time frequency information, and config details
%                    % saved for later retrieval.

%% Get information
sampledur = (channels.time{1}(2)-channels.time{1}(1)); % Duration of each sample
numchans = size(channels.trial{1},1); % Number of carriers
minfft = config.fft(1);
maxfft = config.fft(2);
minfoi = config.foi(1);
maxfoi = config.foi(2);
mintime = config.time(1);
maxtime = config.time(2);
width = config.waveletwidth;

if strcmp(config.cutmethod,'consistenttime')
    mintime_fft = mintime;
    maxtime_fft = maxtime;
elseif strcmp(config.cutmethod,'cutartefact') % cut additional time
    mintime_fft = mintime-0.5; 
    maxtime_fft = maxtime+0.5;
end

%% Calculate FFT
cfgtf           = [];
cfgtf.method    = 'wavelet';
cfgtf.width     = width;
cfgtf.toi       = mintime_fft:sampledur:maxtime_fft;
cfgtf.foi       = (minfft:1:maxfft);
cfgtf.output    = 'fourier';
fspec           = ft_freqanalysis(cfgtf,channels);
fspec.powspctrm = abs(fspec.fourierspctrm);

% Apply temporary 1/F subtraction
fspec_old       = fspec;
fspec           = uni_subtract1f(fspec); % apply 1/F subtraction. This is just temporary to find the right channel

%% Find phase of chans in freq range of interest
minfoi_ind = find(abs(minfoi-fspec.freq)==min(abs(minfoi-fspec.freq))); % minimun freq of interest index
maxfoi_ind = find(abs(maxfoi-fspec.freq)==min(abs(maxfoi-fspec.freq))); % maximum freq of interest index
foivec_ind = minfoi_ind:maxfoi_ind; % freq range of interest vector
mintime_ind = find(abs(mintime-fspec.time)==min(abs(mintime-fspec.time))); % minimun time of interest index
maxtime_ind = find(abs(maxtime-fspec.time)==min(abs(maxtime-fspec.time))); % maximun time of interest index

for chan = 1:numchans
    powtf(:,:,chan)=squeeze(nanmean(fspec.powspctrm(:,chan,:,:),1));
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
    topo = zeros(numel(channels.topolabel));
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
if isfield(config,'topchan')
    numtopchans = min(config.topcarrier,numchans);
else
    numtopchans = min(numchans,30);
end
chanrank = chanrank(1:numtopchans,:); % take the top numcarrier carriers

%% Save information
% Filter relevant frequency spectrum information
fspecinfo.freq = fspec.freq;
fspecinfo.time = fspec.time;

fft_chans{1} = chanrank; %information about top channels
fft_chans{2} = [mintime_ind maxtime_ind];
fft_chans{3} = [minfft maxfft];
fft_chans{4} = fspecinfo; %fft information about channels
fft_chans{5} = powtf(:,:,chanrank(:,1));
fft_chans{6} = pspec(:,chanrank(:,1));
fft_chans{7} = phs(chanrank(:,1)); %phase of the top channels
fft_chans{8} = config.cutmethod;
end
