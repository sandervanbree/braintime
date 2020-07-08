function [fft_comp] = bt_analyzecomps(config, comp)
% Analyze the time frequency characteristics of all ICA components and
% sort them based on average power in the frequency range of interest.
% Optionally, this sorting can be biased by components' correlation
% to a template topography created using bt_temptopo.
%
% Use:
% [fft_comp] = bt_analyzecomps(cfg,comp)

%% Get basic information
sampledur = (comp.time{1}(2)-comp.time{1}(1)); % Duration of each sample
numcomp = size(comp.trial{1},1); % Number of components
minfft = config.minfft;
maxfft = config.maxfft;
minfoi = config.minfoi;
maxfoi = config.maxfoi;
mintime = config.mintime;
maxtime = config.maxtime;

if strcmp(config.cutmethod,'consistenttime')
    mintime_fft = mintime;
    maxtime_fft = maxtime;
elseif strcmp(config.cutmethod,'cutartefact') % Add an additional second that will later be cut
    mintime_fft = config.mintime-0.5;
    maxtime_fft = config.maxtime+0.5;
end

%% Calculate FFT
cfgtf           = [];
cfgtf.method    = 'wavelet';
cfgtf.width     = 5;
cfgtf.toi       = mintime_fft:sampledur:maxtime_fft; %Add 0.5s of data at each end, to be cut out later
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
    oscmaxfreq(cmp) = (fspec.freq(foivec_ind(maxpowind))); % what´s the highest power oscillation?
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

% Filter relevant frequency spectrum information
fspecinfo.freq = fspec.freq;
fspecinfo.time = fspec.time;

%save basic information
fft_comp{1} = comprank; %information about top components
fft_comp{2} = [mintime_ind maxtime_ind];
fft_comp{3} = [minfft maxfft];
fft_comp{4} = fspecinfo; %fft information about components
fft_comp{5} = powtf(:,:,comprank(:,1));
fft_comp{6} = pspec(:,comprank(:,1));
fft_comp{7} = phs(comprank(:,1)); %phase of the top components
fft_comp{8} = config.cutmethod;
end
