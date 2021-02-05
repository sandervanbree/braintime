function [fft_sources] = bt_analyzesources(config,configFT,warpsources)
% Analyze the time frequency characteristics of warping sources (e.g.
% ICA components, virtual channels, LFP time series). The warping sources
% will be ranked according to power at frequencies of interest, and
% may optionally be biased by their topography (see config.rankmethod).
% Each warping source contains a candidate warping signal (the oscillation
% with the most power in a range of interest).
%
% Use:
% [fft_sources] = analyzesources(cfg,cfg.FT,warpsources)
%
% Input Arguments:
% config             % TOOLBOX configuration structure with cells:
%                    %
%   - time           % Start and end of the time window of interest
%   - warpfreqs      % Lowest and highest frequency of interest for the to
%                    % be designated warping signal
%   - correct1f      % (Optional). 'yes' to apply 1/F correction to
%                    % visualized power spectrum (default: 'yes').
%   - nwarpsources   % Number of best warping sources to be extracted from
%                    % the total amount of available warping sources.
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
%   - rankmethod     % 'maxpow': rank warping sources according to their
%                    % average power between minfoi and maxfoi.
%                    %
%                    % 'templatetopo': loads a template topography created
%                    % using bt_templatetopo that is used to bias the
%                    % power ranking to the each source's topography
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
% warpsources        % FieldTrip data structure that contains warping
%                    % sources (ICA components, virtual channels,
%                    % or LFP). Each warping source contains a 
%                    % warping signal, one of which will be selected later.
%                    %
% Output:            %
% fft_sources        % Data structure with: ranked warping sources, their 
%                    % time frequency information, and config details


%% Get information
minfoi = config.warpfreqs(1);                           % Lowest freq of interest
maxfoi = config.warpfreqs(2);                           % Highest freq of interest
mintime = config.time(1);                               % Minimum time of interest
maxtime = config.time(2);                               % Maximum time of interest
cfgFT = configFT;                                       % FieldTrip config structure
correct1f = config.correct1f;                           % Correct for 1/f in visualization of FFT option
nwarpsources = config.nwarpsources;                     % Number of top sources to be considered

try % For ICA components, this should work
sampledur = (warpsources.time{1}(2)-warpsources.time{1}(1));  % Duration of each sample
numsrc = size(warpsources.trial{1},1);                  % Number of warping sources
catch % For virtual channels, this should work
sampledur = (warpsources.time(2)-warpsources.time(1));  % Duration of each sample
numsrc = size(warpsources.trial,2);                     % Number of warping sources  
end

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
fspec           = ft_freqanalysis(cfgFT,warpsources);

fspec.powspctrm = abs(fspec.fourierspctrm);
fspec_old       = fspec;

if strcmp(correct1f,'yes')
fspec           = uni_subtract1f(fspec); % apply 1/F subtraction. This is just temporary to find the right warping source
end

%% Find phase of warping sources in freq range of interest
% Find indices of interest
minfoi_ind = nearest(fspec.freq,minfoi); % index of lowest warping frequency
maxfoi_ind = nearest(fspec.freq,maxfoi); % index of highest warping frequency
foivec_ind = minfoi_ind:maxfoi_ind;      % freq range of interest vector
mintime_ind = nearest(fspec.time,mintime); % minimum time of interest index
maxtime_ind = nearest(fspec.time,maxtime); % maximum time of interest index

% Pre-allocate
powtf = zeros(size(fspec.powspctrm,3),maxtime_ind+1-mintime_ind,numsrc);
pspec = zeros(size(fspec.powspctrm,3),numsrc);
oscmaxfreq = zeros(numsrc,1);
foi_ind = zeros(numsrc,1);
phs = cell(numsrc,1);

% For each source, obtain time frequency information
for src = 1:numsrc
    powtf(:,:,src)=squeeze(nanmean(fspec.powspctrm(:,src,:,mintime_ind:maxtime_ind),1));
    pspec(:,src)=squeeze(nanmean(powtf(:,:,src),2)); % get power spectrum of warping sources
    [maxpow, maxpowind]= max(pspec(minfoi_ind:maxfoi_ind,src)); % what's the highest power in the freq range of interest and its index?
    oscmaxfreq(src) = (fspec.freq(foivec_ind(maxpowind))); % what's the warping signal?
    [~, freqidx]=find(abs(oscmaxfreq(src)-fspec.freq)==min(abs(oscmaxfreq(src)-fspec.freq))); % what is the index of the highest power oscillation?
    foi_ind(src) = freqidx;
    oscmax(src,:) = [src,oscmaxfreq(src), foi_ind(src),maxpow]; % matrix with warping source number, freq, its index, and its power
    phs{src}=squeeze(angle(fspec_old.fourierspctrm(:,src,foi_ind(src),:))); % fetch phase of warping signal
end

%% rank warping sources by one of two methods
% For maxpow this is sufficient:
srcrank = sortrows(oscmax,4,'descend'); % rank sources from highest to lowest power

% Template topography builds onto the variable srcrank
if strcmp(config.rankmethod,'templatetopo')
    try
        load temptopo
        disp('########################################################################################');
        disp('A template topography detected on the toolbox path. If this is from a previous analysis,')
        disp('please rerun bt_templatetopo and create a new template topography');
        disp('########################################################################################');
    catch
        error('A template topography was not found. Use bt_templatetopo to create a template topography')
    end
    
    % create final template topography
    try
        topo_ind = contains(warpsources.topolabel,temptopo); % Find the indices of the chosen sensors in the channels' labels
    catch
        error('Channel labels in your template topography mismatch carrier labels')
    end
    topo = double(topo_ind)*5; % Put high value (5) only for chosen channels
    
    % Get the source-template correlation
    corrord=srcrank;
    corrord(:,1)=1:numsrc;
    
    % correlate template topography with source topography
    for src=1:numsrc
        tc=corrcoef(topo,abs(warpsources.topo(:,src)));
        corrord(src,2)=tc(1,2);
    end
    corrord=sortrows(corrord,2,'descend'); % rank sources from highest to lowest correlation
    corrord=corrord(:,1);
    
    crank=numsrc;
    for i=1:crank
        srcrank(i,5)=crank; % add descending values to the power order sources
        corrord(i,2)=crank; % add descending values to the correlation order sources
        crank=crank-1;
    end
    srcrank=sortrows(srcrank,1,'ascend'); % rank by source number
    corrord=sortrows(corrord,1,'ascend');
    
    srcrank(:,6)=srcrank(:,5)+corrord(:,2); %get an index with the addition
    
    srcrank=sortrows(srcrank,6,'descend'); % rank sources by index
    srcrank=srcrank(:,1:4);
end
    
%% Take the best sources (default: 30)
if isfield(config,'nwarpsources')
    nwarpsources = min(nwarpsources,numsrc);
else
    nwarpsources = min(numsrc,30);
end
srcrank = srcrank(1:nwarpsources,:); % take the best sources

%% Save information
% Filter relevant frequency spectrum information
fspecinfo.freq = fspec.freq;
fspecinfo.time = fspec.time;

fft_sources{1} = srcrank;                                  % Time freq data of extracted warping sources
fft_sources{2} = [mintime maxtime];                        % Start and end time of interest
fft_sources{3} = [minfoi maxfoi];                          % Lowest and highest possible warping frequency
fft_sources{4} = fspecinfo;                                % FFT time and frequency vector
fft_sources{5} = powtf(:,:,srcrank(:,1));                  % Power spectrum of warping sources
fft_sources{6} = pspec(:,srcrank(:,1));                    % Power spectrum averaged across trials
fft_sources{7} = phs(srcrank(:,1));                        % Phase of all warping sources for all trials
fft_sources{8} = config.cutmethod;                         % Applied cutting method
fft_sources{9} = config.rankmethod;                       % Applied ranking method
end