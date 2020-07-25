function [bt_carrier] = bt_GEDanalyzechoose(config, channels)
% Documentation to follow
% Adapted from Mike Cohen, 2017; DOI: 10.7554/eLife.21792

%% Get information
sampledur = (channels.time{1}(2)-channels.time{1}(1)); % Duration of each sample
samplerate = 1/sampledur;
trialnum = size(channels.trial,2); % number of trials
numchans = size(channels.label,1);
samples = size(channels.trial{1},2);
minfoi = config.foi(1);
maxfoi = config.foi(2);
foi = maxfoi-(maxfoi-minfoi)/2;
mintime = config.time(1);
maxtime = config.time(2);
duration = maxtime-mintime;
Ntopchan = config.Ntopchan;
fwhm = config.fwhm;
cutmethod = config.cutmethod;
width = config.waveletwidth;
removecomp = config.removecomp;

if strcmp(config.cutmethod,'consistenttime')
    mintime_fft = mintime;
    maxtime_fft = maxtime;
    mintime_ind = 1;
    maxtime_ind = samples;
elseif strcmp(config.cutmethod,'cutartefact') % cut additional time
    mintime_fft = mintime-0.5;
    maxtime_fft = maxtime+0.5;
    mintime_ind = round((samplerate*0.5));
    maxtime_ind = round((maxtime_fft-mintime_fft)*samplerate-samplerate*0.5);
end

if ~isfield(config,'fwhm')
    error('For gedcomp, choose a full width half maximum using config.fwhm')
end

% Cut trial channels for peak finding
cfg = [];
cfg.toilim = [mintime maxtime];
channels_pk = ft_redefinetrial(cfg,channels);

%% Find peak
[PSD,f]=pwelch(channels_pk.trial',samplerate,0,[],samplerate,'psd');
PSD=mean(bsxfun(@rdivide,PSD,sum(PSD,1)),2);
PSD=PSD(f>minfoi&f<maxfoi,1); f=f(f>minfoi&f<maxfoi);
[TREND]=polyfit(log10(f),log10(PSD),3);
Y=polyval(TREND,log10(f));
fpeak=find((log10(PSD)-Y)==max(log10(PSD)-Y));
pk=f(fpeak);
freqpk=round(pk,1);
disp(['Continuing with peak in frequency range of interest identified at ' num2str(round(pk,1)) ' Hz.'])

%Cut trials for GED comp
cfg = [];
cfg.toilim = [mintime_fft maxtime_fft];
channels_GED = ft_redefinetrial(cfg,channels);

%% GED to identify theta component
osccov   = zeros( size(channels_GED.label,1) ); % preallocate covariance matrix of filtered data
bbcov    = zeros( size(channels_GED.label,1) ); % preallocate covariance matrix of broadband data
trialnum = size(channels_GED.trial,2); % number of trials
oscfilt  = zeros(trialnum, size(channels_GED.label,1), size(channels_GED.trial{1},2));

for trl = 1: trialnum % loop through trials    
    oscfiltTrl = filterFGx(channels_GED.trial{trl},channels_GED.fsample,freqpk,fwhm); % filter data in frequency of interest
    oscfiltTrl = bsxfun(@minus,oscfiltTrl,mean(oscfiltTrl,2)); % demean (this is important because otherwise the first component will likely point towards the mean of the covariance matrix)
    tempcov = (oscfiltTrl*oscfiltTrl')/samples; % covariance matrix of the filtered signal for trial trl; the covariance is between channels
    osccov  = osccov + tempcov; % add this to the total covariance matrix over all trials
    oscfilt(trl,:,:) = oscfiltTrl;
    
    % broadband covariance
    tmpdat = bsxfun(@minus,channels_GED.trial{trl},mean(channels_GED.trial{trl},2)); % demean broadband data
    tempcov  = (tmpdat*tmpdat')/samples; % compute broadband covariance matrix
    bbcov = bbcov+tempcov;    
end

% divide by number of trials (this might not be necessary, I'm not sure
% right now)
osccov = osccov / trialnum;
bbcov = bbcov / trialnum;

% GED
% compute eigenvalues between alpha covariance matrix and broadband covariance matrix; this is not a normal eigenvaluedecomposition, but a generalized form.
% you find it very well explained in the mike cohen paper.
% the biggest eigenvector points towards the highest quotient of the two input covariance matrices.
[evecsT,evals] = eig(osccov,bbcov);

[~,chanrank] = sort(diag(evals),'descend'); % index of the biggest component; usually it does not need sorting, but apparently it is possible
chanrank = chanrank(1:Ntopchan);

% multiply the biggest eigenvector with the transpose of the filtered channels.
% you can understand this as a weighted average of the channels.
% You can also get a topography (I think by multiplying the eigenvector with the covariance matrix (but see method1 in the paper to be sure)
mat     = permute(oscfilt,[2 3 1]);
allComp = reshape(mat, size(oscfilt,2), [] );

uiwait(msgbox({'The components are sorted based on the presence of the identified frequency of interest peak';' ';...
    'Please pick one component, factoring in its:';...
    '(1) Time series; when do you expect high power?';...
    '(2) Component number; the lower the component number the higher the explained r^2.';' ';...
    'Note: topography plotting is not currently implemented for this method, so use bt_analyzechannels and bt_choosecarrier instead if you want to choose a carrier based on regions of interest'}))
uiwait(msgbox({'Instructions to browse through components:';' ';...
    'Press forward/back arrow to see the next/previous component';...
    'Once you have decided for one component, click on that component and press Q to quit visualization'}))

chanind = 1;
finish = 0;
while finish==0
    if chanind < 1 % make sure channel cannot go out of bounds
        chanind = 1;
    elseif chanind>size(chanrank,1)
        chanind = size(chanrank,1);
    end
    
    osccomp   = allComp' * evecsT(:,chanrank(chanind));
%     if strcmp(cutmethod,'cutartefact')
%         osccomp = osccomp(trialnum*startind:trialnum*endind);
%     end
    
    figure(1)
    clf
    
    % component
    subplot(211)
    plot(osccomp, 'linew', 1)
    xlim([1 numel(osccomp)])
    xticks(linspace(1,numel(osccomp),10));
    xticklabels(round(linspace(mintime,maxtime,10),2));
    ylim([-max(abs(get(gca,'YLim'))) max(abs(get(gca,'YLim')))])
    xlabel('Time (s)')
    ylabel('Component Amplitude')
    title(sprintf('GED Component Time Series | Component %d at peak %0.2fHz',chanind,freqpk))
    
    
    % powerspectrum
    subplot(224)
    for trl = 1:trialnum
        temp_PS = squeeze(oscfilt(trl,:,mintime_ind:maxtime_ind))' * evecsT(:,chanrank(chanind));
        temp_PHS = squeeze(oscfilt(trl,:,:))' * evecsT(:,chanrank(chanind));
        powSpec(trl,:) = abs(fft(temp_PS)); %The power spectrum is for plotting, so restrict to toi
        phs(trl,:) = angle(fft(temp_PHS)); %The phase needs to be longer for cutartefact
    end
    
    powSpec = mean(powSpec,1);
    
    hz = linspace(0, channels.fsample, size(powSpec,2));
    plot(hz,powSpec, 'linew', 2)
    xlim([0 30])
    xlabel('Frequency')
    ylabel('Power')
    title(' GED Component Powerspectrum')
    
    
    keydown = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    [~,~,keyCode] = KbCheck;
    key = KbName(find(keyCode));
    if (keydown == 0) %grab the component if click
        channeloi=chanind;
        fprintf('Selected component number %d. Press ''q'' to quit.\n',chanind)
        chosenVec = evecsT(:,chanrank(chanind));
    elseif value == 28 % left arrow
        if chanind == 1
            disp('This is the first component')
        else
            chanind = chanind-1;
        end
        fprintf('Component #%d\n', chanind);
    elseif value == 29 % right arrow
        if chanind == numel(chanrank) % make sure channel cannot go out of bounds
            disp('This is the last component')
        else chanind = chanind + 1;
        end
        fprintf('Component #%d\n', chanind);
    elseif strcmp(key,'q') %stop the loop if it is not necessesary to keep visualising
        fprintf('Carrier will be the phase in component number %d.\n',chanind)
        chanind = (numel(chanrank))+1;
        finish = 1;
        close
    end
end


for trl = 1:trialnum
    oscMC            =  squeeze(oscfilt(trl,:,:))   .* chosenVec; % component per trial (1xnpnts)
    oscComp(trl,:)   =  squeeze(oscfilt(trl,:,:))'   * chosenVec; % each row is a trial each column is a timpoint within eahc trial
    
    if isfield(config,'removecomp')
        if strcmp(config.removecomp,'no')
            break
        else
%             channels.trial{trl} =  channels.trial{trl} - oscMC; % remove component from original data
        end        
    end    
end



% Filter relevant frequency spectrum information
fspecinfo.freq = hz;
fspecinfo.time = linspace(mintime_fft,maxtime_fft,numel(hz));

bt_carrier{1} = channeloi; %chosen carrier
bt_carrier{2} = {phs}; %phase of chosen carrier
bt_carrier{3} = channels;
bt_carrier{4} = [channeloi,freqpk];
bt_carrier{5} = fspecinfo;
bt_carrier{6} = cutmethod;
bt_carrier{7} = [mintime_ind, maxtime_ind];
bt_carrier{8} = 'bt_GEDanalyzechoose'; %Label method used
end

