function [fft_comp] = bt_analyzecomps(config, comp)
% Analyze the time frequency characteristics of all ICA components and
% sort them based on average power in the frequency range of interest.
% Optionally, this sorting can be biased by components' correlation
% to a template topography created using bt_temptopo.
%
% Use:
% [fft_comp] = bt_analyzecomps(cfg,comp)
%
% Input Arguments:
% config
%   - mintime        % Start time window of interest
%   - maxtime        % End time window of interest
%   - minfft         % Lowest frequency analyzed in the components 
%   - maxfft         % Highest frequency analyzed in the components
%   - minfoi         % Lowest frequency of interest for the to be
%                    % designated brain time oscillation (e.g. low alpha
%                    % for attention)
%   - maxfoi         % Highest frequency of interest for the to be
%                    % designated brain time oscillation (e.g. high alpha
%                    % for attention)
%   - topcomp        % Number of best components to be filtered from the full amount
%                    
%   - cutmethod      % 'consistenttime': warp from mintime to maxtime.
%                    % The upside of this method is that the final
%                    % brain time data is equally long and of the same
%                    % duration as intended (mintime to maxtime). The
%                    % downside is an artefact in the first cycle
%                    % caused by the warping requiring a repetition of data
%                    % bins for alignment to the template oscillation.
%                    % 
%                    % 'cutartefact': warp from mintime-0.5sec to
%                    % maxtime+0.5sec, later cut to mintime maxtime.
%                    % The upside of this method is that the final brain
%                    % time data has no artefact in the first cycle. The
%                    % downside is that there is variance across trials
%                    % in the brain time's start and end time, as well as
%                    % the duration.
%                             
%   - sortmethod     % 'maxpow': sort components according to their average
%                    % power in minfoi and maxfoi.
%                    %
%                    % 'temptopo': loads a template topography created
%                    % using bt_templatetopo that is used to bias the
%                    % power sorting to the each component's topography
%                    % match to the template topography.
%                    %
%   - removecomp     % 'yes': removes component from the brain time data.
%                    % When analyzing brain time data using your own
%                    % analysis pipeline, you may wish to remove the 
%                    % component to avoid circularity. See the brain time
%                    % toolbox paper for more details.
%                    %
%                    % 'no': keeps component in the brain time data.
%                    %
% comp               % FieldTrip component data structure as obtained
%                    % by applying ft_componentanalysis on clock time data.
%                    %
% Output:            %
% fft_comp           % Data structure with: ranked components, their 
%                    % time frequency information, and config details
%                    % saved for later retrival.

%% Get information
sampledur = (comp.time{1}(2)-comp.time{1}(1)); % Duration of each sample
numcomp = size(comp.trial{1},1); % Number of components
minfft = config.minfft;
maxfft = config.maxfft;
minfoi = config.minfoi;
maxfoi = config.maxfoi;
mintime = config.mintime;
maxtime = config.maxtime;
extracut = 1/minfoi; %cut slowest freq + two samples

if strcmp(config.cutmethod,'consistenttime')
    mintime_fft = mintime;
    maxtime_fft = maxtime;
elseif strcmp(config.cutmethod,'cutartefact') % cut additional time
    mintime_fft = config.mintime-extracut; 
    maxtime_fft = config.maxtime+extracut;
end

%% Calculate FFT
cfgtf           = [];
cfgtf.method    = 'wavelet';
cfgtf.width     = 5;
cfgtf.toi       = mintime_fft:sampledur:maxtime_fft;
cfgtf.foi       = (minfft:1:maxfft);
cfgtf.output    = 'fourier';
fspec           = ft_freqanalysis(cfgtf,comp);
fspec.powspctrm = abs(fspec.fourierspctrm);

% Apply temporary 1/F subtraction
fspec_old       = fspec;
fspec           = uni_subtract1f(fspec); % apply 1/F subtraction. This is just temporary to find the right component

%% Find phase of components in freq range of interest
minfoi_ind = find(abs(minfoi-fspec.freq)==min(abs(minfoi-fspec.freq))); % minimun freq of interest index
maxfoi_ind = find(abs(maxfoi-fspec.freq)==min(abs(maxfoi-fspec.freq))); % maximum freq of interest index
foivec_ind = minfoi_ind:maxfoi_ind; % freq range of interest vector
mintime_ind = find(abs(mintime-fspec.time)==min(abs(mintime-fspec.time))); % minimun time of interest index
maxtime_ind = find(abs(maxtime-fspec.time)==min(abs(maxtime-fspec.time))); % maximun time of interest index

for cmp = 1:numcomp
    powtf(:,:,cmp)=squeeze(nanmean(fspec.powspctrm(:,cmp,:,mintime_ind:maxtime_ind),1));
    pspec(:,cmp)=squeeze(nanmean(powtf(:,:,cmp),2)); % get power spectrum of components
    [maxpow, maxpowind]= max(pspec(minfoi_ind:maxfoi_ind,cmp)); % what's the highest power in the freq range of interest and its index in the freq range of interest vector?
    oscmaxfreq(cmp) = (fspec.freq(foivec_ind(maxpowind))); % what�s the highest power oscillation?
    [~, freqidx]=find(abs(oscmaxfreq(cmp)-fspec.freq)==min(abs(oscmaxfreq(cmp)-fspec.freq))); % what is the index of the highest power oscillation?
    foi_ind(cmp) = freqidx;
    oscmax(cmp,:) = [cmp,oscmaxfreq(cmp), foi_ind(cmp),maxpow]; % matrix with component number, freq, its index, and its power
    phs{cmp}=squeeze(angle(fspec_old.fourierspctrm(:,cmp,foi_ind(cmp),:))); % fetch phase of highest power oscillation
end

%% Sort components by one of two methods
% For maxpow this is sufficient:
comprank=sortrows(oscmax,4,'descend'); % sort components from highest to lowest power

% Template topography builds onto the variable comprank
if strcmp(config.sortmethod,'temptopo')
    try
        load temptopo
    catch
        error('Template topography not found. Use bt_templatetopo to create a template topography')
    end
    
    % create final template topography
    topo = zeros(numel(comp.topolabel));
    try
        topo_ind = contains(comp.topolabel,temptopo); % Find the indices of the chosen channels in the component's labels
    catch
        error('Channel labels in your template topography mismatch component labels')
    end
    topo = double(topo_ind)*5; %Put high value (5) only for chosen channels
    
    % Get the components-template correlation
    corrord=comprank;
    corrord(:,1)=1:numcomp;
    
    % correlate template topography with component topography
    for cmp=1:numcomp
        tc=corrcoef(topo,abs(comp.topo(:,cmp)));
        corrord(cmp,2)=tc(1,2);
    end
    corrord=sortrows(corrord,2,'descend'); % sort components from highest to lowest correlation
    corrord=corrord(:,1);
    
    crank=numcomp;
    for i=1:crank
        comprank(i,5)=crank; % add descending values to the power order components
        corrord(i,2)=crank; % add descending values to the correlation order components
        crank=crank-1;
    end
    comprank=sortrows(comprank,1,'ascend'); % sort by components number
    corrord=sortrows(corrord,1,'ascend');
    
    comprank(:,6)=comprank(:,5)+corrord(:,2); %get an index with the addition
    
    comprank=sortrows(comprank,6,'descend'); % sort components by index
    comprank =comprank(:,1:4);
end

%% Take the 30 best components (or specified number)
if isfield(config,'topcomps')
    numtopcomps = min(config.topcomps,numcomp);
else
    numtopcomps = min(numcomp,30);
end
comprank = comprank(1:numtopcomps,:); % take the top numcomp components

%% Save information
% Filter relevant frequency spectrum information
fspecinfo.freq = fspec.freq;
fspecinfo.time = fspec.time;

fft_comp{1} = comprank; %information about top components
fft_comp{2} = [mintime_ind maxtime_ind];
fft_comp{3} = [minfft maxfft];
fft_comp{4} = fspecinfo; %fft information about components
fft_comp{5} = powtf(:,:,comprank(:,1));
fft_comp{6} = pspec(:,comprank(:,1));
fft_comp{7} = phs(comprank(:,1)); %phase of the top components
fft_comp{8} = config.cutmethod;
fft_comp{9} = extracut;
end
