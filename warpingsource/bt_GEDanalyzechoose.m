function [bt_carrier, channels] = bt_GEDanalyzechoose(config, channels)
% Documentation to follow
% Adapted from Mike Cohen, 2017; DOI: 10.7554/eLife.21792
% Not implemented yet: topography plotting, component removal

% DEV NOTES ~ CURRENT ISSUES I NEED HELP WITH:
% ISSUE 1:
% Unlike ICA, warping the data with the phase of the best component as obtained by GED
% does not give recurrence when you run the obtained bt_carrier structure
% through tutorial 2. This most likely means I have done something wrong
% in my implementation of this method. To find the mistake, I can recommend
% ignoring the method 'cutartefact' because it may cause time window
% confusion. 'consistenttime' is without any change in time at any point.
%
% ISSUE 2:
% Mike Cohen's tools to plot the GED components' resources requires an 
% EEG Lab layout structure that we are not using or familiar with.
% I am unsure how to adapt that script ("topoplotIndie"); see https://github.com/cogneurd/GED
% for it to work with the toolbox

% ISSUE 3: 
% I am unsure how to do component removal because to obtain component time
% courses, we have to cut to our time window of interest (otherwise the
% component optimization will base itself of time frequency info outside
% the window of interest). But that also means that the variable 
% we would ordinarily subtract from the data ("oscMC") is of the wrong
% data sample length.

%% Get information
sampledur = (channels.time{1}(2)-channels.time{1}(1));    % Duration of each sample
sr = 1/sampledur;                                         % Sampling rate
trialnum = size(channels.trial,2);                        % Number of trials
minfoi = config.foi(1);                                   % Lowest freq of interest
maxfoi = config.foi(2);                                   % Highest freq of interest
mintime = config.time(1);                                 % Minimum time of interest
maxtime = config.time(2);                                 % Maximum time of interest
Ntopchan = config.Ntopchan;                               % Number of best components to be considered
fwhm = config.fwhm;                                       % Full width at half maximum to be used for peak picking
cutmethod = config.cutmethod;                             % Applied cutting method

% Adapt time for the two cutting methods
if strcmp(config.cutmethod,'consistenttime')
    mintime_fft = mintime;
    maxtime_fft = maxtime;
elseif strcmp(config.cutmethod,'cutartefact') % cut additional time
    mintime_fft = mintime-0.5;
    maxtime_fft = maxtime+0.5;
end

if ~isfield(config,'fwhm')
    error('For gedcomp, choose a full width half maximum using config.fwhm')
end

% Cut trial channels for peak finding
cfg = [];
cfg.toilim = [mintime maxtime];
channels_pk = ft_redefinetrial(cfg,channels);

%% Find peak in channel data
[PSD,f]=pwelch(channels_pk.trial',sr,0,[],sr,'psd'); % Use p-welch method to find power spectrum
PSD=mean(bsxfun(@rdivide,PSD,sum(PSD,1)),2);
PSD=PSD(f>minfoi&f<maxfoi,1); f=f(f>minfoi&f<maxfoi); % Filter in freq range of interest
[TREND]=polyfit(log10(f),log10(PSD),3); % Fit 1/f pattern
Y=polyval(TREND,log10(f));
fpeak=find((log10(PSD)-Y)==max(log10(PSD)-Y)); % Find frequency peak
pk=f(fpeak);
freqpk=round(pk,1); %This is the frequency that will be used for GED
disp(['Continuing with peak in frequency range of interest identified at ' num2str(round(pk,1)) ' Hz.'])

%Cut trials for GED comp
cfg = [];
cfg.toilim = [mintime_fft maxtime_fft];
channels_GED = ft_redefinetrial(cfg,channels);
mintime_ind = find(abs(mintime-channels_GED.time{1})==min(abs(mintime-channels_GED.time{1}))); % minimun time of interest index
maxtime_ind = find(abs(maxtime-channels_GED.time{1})==min(abs(maxtime-channels_GED.time{1}))); % maximun time of interest index
samples = size(channels_GED.trial{1},2);

%% GED to identify component of frequency of interest
% Pre-allocate
osccov   = zeros( size(channels_GED.label,1) ); 
bbcov    = zeros( size(channels_GED.label,1) );
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

% Divide by number of trials
osccov = osccov / trialnum;
bbcov = bbcov / trialnum;

% GED: compute eigenvalues between alpha covariance matrix and broadband covariance matrix
[evecsT,evals] = eig(osccov,bbcov); % The biggest eigenvector points towards the highest quotient of the two input covariance matrices.

% Sort from high to low
[~,chanrank] = sort(diag(evals),'descend'); %Note: In contrast to ICA, higher component numbers have a higher weight
chanrank = chanrank(1:Ntopchan); %Filter on top channels

% Reshape data into workable format
mat     = permute(oscfilt,[2 3 1]);
allComp = reshape(mat, size(oscfilt,2), [] );

% % Pre-allocate
powSpec = zeros(1,numel(mintime_ind:maxtime_ind));
phs = zeros(trialnum,numel(mintime_ind:maxtime_ind));
oscMC = zeros(size(oscfilt,2),size(oscfilt,3));
oscComp = zeros(size(oscfilt,3));

uiwait(msgbox({'The components are sorted based on the presence of the identified frequency of interest peak';' ';...
    'Please pick one component, factoring in its:';...
    '(1) Time series; when do you expect high power?';...
    '(2) Component number; the lower the component number the higher the component''s weight.';' ';...
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
    if strcmp(cutmethod,'cutartefact')
        osccomp = osccomp(trialnum*mintime_ind:trialnum*maxtime_ind);
    end
    
    figure(1)
    clf
    
    % Plot component time series
    subplot(211)
    plot(osccomp, 'linew', 1)
    xlim([1 numel(osccomp)])
    xticks(linspace(1,numel(osccomp),10));
    xticklabels(round(linspace(mintime,maxtime,10),2));
    ylim([-max(abs(get(gca,'YLim'))) max(abs(get(gca,'YLim')))])
    xlabel('Time (s)')
    ylabel('Component Amplitude')
    title(sprintf('GED Component Time Series | Component %d at peak %0.2fHz',chanind,freqpk))
    
    % Plot the power spectrum
    subplot(224)
    for trl = 1:trialnum
        temp_PS = squeeze(oscfilt(trl,:,mintime_ind:maxtime_ind))' * evecsT(:,chanrank(chanind)); % Only take time window of interest
        temp_PHS = squeeze(oscfilt(trl,:,:))' * evecsT(:,chanrank(chanind)); % But get the full time for phase
        powSpec(trl,:) = abs(fft(temp_PS)); % The power spectrum is for plotting
        phs(trl,:) = angle(fft(temp_PHS)); % Get the phase
    end
    
    % Average the power spectrum
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
        else
            chanind = chanind + 1;
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
    oscComp(trl,:)   =  squeeze(oscfilt(trl,:,:))'   * chosenVec; % each row is a trial each column is a timepoint within each trial
    
    if isfield(config,'removecomp')
        if strcmp(config.removecomp,'no')
            break
        else
%          channels.trial{trl} =  channels.trial{trl} - oscMC; % remove component from original data
        error('Component removal not yet implemented in GED method')
        end        
    end    
end

% Filter relevant frequency spectrum information
fspecinfo.freq = hz;
fspecinfo.time = linspace(mintime_fft,maxtime_fft,numel(hz));

bt_carrier{1} = channeloi;                          % Time freq data of chosen component
bt_carrier{2} = {phs};                              % Phase of chosen component for all trials
bt_carrier{3} = channels;                           % Channel data
bt_carrier{4} = [channeloi,freqpk];                 % Time freq data of chosen channel
bt_carrier{5} = fspecinfo;                          % FFT time and frequency vector
bt_carrier{6} = cutmethod;                          % Applied cutting method
bt_carrier{7} = [mintime_ind, maxtime_ind];         % Index of start and end time of interest (differs for cutartefact)
bt_carrier{8} = 'bt_GEDanalyzechoose';              % Save carrier choosing method
end

