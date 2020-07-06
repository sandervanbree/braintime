function [fft_comp] = bt_analyzecomps(config, comp)
% Help function to be added here

% Get basic information
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

% Calculate FFT
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

% Find phase of components in freq range of interest
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

% Take the 30 best components (or specified number)
if isfield(config,'topcomps')
numtopcomps = min(config.topcomps,numcomp);
else
numtopcomps = min(numcomp,30);
end

% Sort components
if strcmp(config.sortmethod,'maxpow')
pspec_comp=sortrows(oscmax,4,'descend'); % sort components from highest to lowest power
topcomps = pspec_comp(1:numtopcomps,:); % take the top ncmp components
elseif strcmp(config.sortmethod,'temptopo') %currently not working yet
load topotemp
temp_corr=foicomps;
for cmp=1:ncmp
   tc=corrcoef(topotemp,abs(comp.topo(:,foicomps(cmp))));
   temp_corr(cmp,2)=tc(1,2);
end
temp_corr=sortrows(temp_corr,2,'descend');
end

% Filter relevant frequency spectrum information
fspecinfo.freq = fspec.freq;
fspecinfo.time = fspec.time;

fft_comp{1} = topcomps; %basic information about top components
fft_comp{2} = [mintime_ind maxtime_ind]; 
fft_comp{3} = [minfft maxfft];
fft_comp{4} = fspecinfo; %basic fft information about components
fft_comp{5} = powtf(:,:,topcomps(:,1));
fft_comp{6} = pspec(:,topcomps(:,1));
fft_comp{7} = phs(topcomps(:,1)); %phase of the top components
fft_comp{8} = config.cutmethod;
end
