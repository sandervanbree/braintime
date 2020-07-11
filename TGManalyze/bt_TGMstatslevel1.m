function [stats1] = bt_TGMstatslevel1(config, bt_data, bt_TGMquant)
% Acquire single subject level statistics using permutation testing.
% Null distributions are created by shuffling the classification labels
% numperms1 times and collecting the power spectra from the resulting AC
% maps. If config.figure = 'yes', displays stats results on the single
% subject level.
% To enable second level testing, create one output structure per
% participant (e.g. TGMstat1{currentsubject}).
% 
% Use:
% [stats1{subj}] = bt_TGMstatslevel1(config, bt_data, bt_TGMquant)
%
% Input Arguments:
% config
%   - mvpacfg        % Load the same cfg file as used to generate empirical
%                    % TGM. Critical: difference in shuffle and empirical
%                    % data should not be caused by differences in
%                    % classification parameters.
%                    %
%   - numperms1      % Number of permutations on the first statistical
%                    % level.
%                    %
%   - statsrange     % Range of recurrence rates to be statistically tested
%                    % in the TGM.
%                    %
%   - figure         % 'yes' or 'no': display statistical results
%                    %
% bt_data            % Data structure obtained from bt_TGMquant. Contains:
%                    % TGM, AC map, FFT information of the AC map, and
%                    % config details saved for later retrieval.
%                    %
% bt_TGMquant        %  TGM obtained by mv_classify_timextime
%                    %
% Output:
% stats1             % Output structure which contains power spectra of
%                    % the empirical data, permutation data, the associated
%                    % frequency vector, and the power at the mode freq.

%% Get information
toi = bt_TGMquant.toi;
warpfreq = bt_TGMquant.warpfreq;
acfft = bt_TGMquant.acfft; % AC map FFT information
modefreq = bt_TGMquant.modefreq;
modefreqind = bt_TGMquant.modefreqind; %Index of the mode frequency
TGM = bt_TGMquant.TGM;
refdimension = bt_TGMquant.refdimension;
timevec = bt_TGMquant.timevec;
normalizer = bt_TGMquant.normalizer;
clabel = config.clabel;
cfg_mv = config.mvpacfg;

if isfield(config,'statsrange') %What is the recurrence rate range of interest?
    statsrange = config.statsrange(1):config.statsrange(2);
else
    statsrange = 1:40;
end

%% statistically test TGM
% Calculate autocorrelation map (AC)
ac=autocorr2d(TGM);

% Size of all rows and columns
nvecs=numel(ac(:,1));

% Perform FFT over one row to get f and find out statsrange indices
[~,f]=Powspek(ac(1,:),nvecs/normalizer);
l = find(abs(statsrange(1)-f)==min(abs(statsrange(1)-f))); %minimum frequency to be tested
h = find(abs(statsrange(end)-f)==min(abs(statsrange(end)-f))); %maximum frequency to be tested
srange = l:h;

% Run FFT over all rows and columns of the AC map
for vec=1:nvecs
    %1st dimenssion
    [PS,f]=Powspek(ac(vec,:),nvecs/normalizer);
    PS1(vec,:) = PS(srange);
    
    %2nd dimension
    [PS,f]=Powspek(ac(:,vec),nvecs/normalizer);
    PS2(vec,:) = PS(srange);
    
end
avg_PS = mean(PS1,1)+mean(PS2,1); %Mean power spectra
modepow_emp = avg_PS(modefreqind); %Mean power at mode freq
fullspec_emp = avg_PS;

% First level permutations
for perm1 = 1:config.numperms1
    fprintf('First level permutation number %i\n', perm1);
    clabel = clabel(randperm(numel(clabel)));
    [permTGM, ~] = mv_classify_timextime(cfg_mv, bt_data.trial, clabel);
    
    % Calculate autocorrelation map (AC)
    ac=autocorr2d(permTGM);
    
    % Run FFT over all rows and columns of the AC map
    nvecs=numel(ac(:,1));
    
    for vec=1:nvecs
        %1st dimenssion
        [PS,f]=Powspek(ac(vec,:),nvecs/normalizer);
        PS1(vec,:) = PS(srange);
        
        %2nd dimension
        [PS,f]=Powspek(ac(:,vec),nvecs/normalizer);
        PS2(vec,:) = PS(srange);
        
    end
    avg_PS = mean(PS1,1)+mean(PS2,1); %Mean power spectra
    modepow_shuff(perm1) = avg_PS(modefreqind); %Mean power spectra at mode freq
    fullspec_shuff(perm1,:) = avg_PS;
end
mean_modepow_shuff = mean(modepow_shuff);

f=f(l:f); %filter frequency vector based on range of interest

% Only calculate confidence interval and plot stats if desired
if isfield(config,'figure')
    if strcmp(config.figure,'yes')
        figopt = 1;
    else
        figopt = 0;
    end
else
    figopt = 1; %Default yes
end

if figopt == 1
    %% Create confidence interval for each frequency bin
    for freq = f
        
        %Grab data
        temp_data = fullspec_shuff(:,freq);
        
        %Calculate mean
        temp_mean = mean(temp_data);
        
        %Calculate standard mean error (SEM)
        temp_SEM = std(temp_data)/sqrt(length(temp_data));
        
        %Calculate t-score
        temp_ts = tinv([0.025  0.975],length(temp_data)-1);
        
        %Calculate confidence interval (CI)
        temp_CI = temp_mean + temp_ts*temp_SEM;
        
        freq_CI(freq,:) = temp_CI; % save CI
    end
    
    %% Plot results
    % 1st plot: relationship empirical and shuffled amp at mode freq
    figure
    subplot(1,2,1)
    x = ones(1,numel(modepow_shuff));
    plot(x,modepow_shuff,'k.','MarkerSize',25,'LineWidth',3);
    grid on;
    hold on;
    plot(x(1),modepow_emp,'r.','MarkerSize',35,'LineWidth',3);
    legend('Shuffled mode power','Empirical mode power')
    ax = gca;
    ax.XTick = 1;
    ax.XTickLabels = {[' ']};
    % Set up axes.
    ylabel('Mean amplitude at mode frequency')
    xlabel(sprintf('%.2f Hz',modefreq))
    xlim([0.8, 1.2]);
    title(sprintf('Amplitude of mode frequency (%.2f Hz)',modefreq))
    
    % 2nd plot: relationship empirical and shuffled amp at all freq
    subplot(1,2,2); hold on
    low_CI = freq_CI(:,1)';
    hi_CI = freq_CI(:,2)';
    p3 = area(f(1,srange),hi_CI);
    p3(1).FaceColor = [0.8627 0.8627 0.8627];
    hold on
    p2 = area(f(1,srange),low_CI);
    p2(1).FaceColor = [1 1 1];
    p1 = plot(f(1,srange),fullspec_emp,'LineWidth',2,'Color','b');
    p2 = plot(f,low_CI,'LineWidth',0.5,'Color','k');
    p3 = plot(f,hi_CI,'LineWidth',0.5,'Color','k');
    
    if strcmp(refdimension,'braintime') %warp freq line is dependent on clock (warped freq) or brain time (1 hz)
        p4 = line([1 1], [0 max(fullspec_emp)],'color',[1 0 1],'LineWidth',4); %Line at warped freq
        xlabel('Recurrence frequency (factor of warped freq)')
    else
        p4 = line([warpfreq warpfreq], [0 max(fullspec_emp)],'color',[1 0 1],'LineWidth',4); %Line at warped freq
        xlabel('Recurrence frequency')
    end
    p4.Color(4) = 0.45;
    
    % set up axes
    ylabel('Mean amplitude')
    title(sprintf('Amplitude across recurrence rates'))
    legend([p1 p2 p3 p4],{'Empirical amplitude','Lower limit CI','Higher limit CI','Warped frequency'});
end

%% Create output structure
stats1.f = f; %frequency spectrum
stats1.empspec = fullspec_emp; %power spectrum of average empirical data
stats1.shuffspec = fullspec_shuff; %power spectrum of average permutation data
stats1.empmode = modepow_emp; %power at mode freq of epirical data
stats1.shuffmode = modepow_shuff; %power at mode freq of permutation data