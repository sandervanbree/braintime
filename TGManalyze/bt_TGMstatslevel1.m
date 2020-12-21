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
%   - normalize      % Normalize empirical and shuffled TGMs by the mean
%                    % and std of shuffled TGMs (Default: yes)
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
numperms1 = config.numperms1;                             % Number of first level permutations
warpfreq = bt_TGMquant.warpfreq;                          % Warped frequency (frequency of the carrier)
modefreq = bt_TGMquant.modefreq;                          % The mode frequency across all rows and columns of the AC map
TGM = bt_TGMquant.TGM;                                    % Time Generalization Matrix of the data
refdimension = bt_TGMquant.refdimension;                  % Reference dimension used
clabel = config.clabel;                                   % Classification labels
cfg_mv = config.mvpacfg;                                  % MVPA Light configuration structured used to obtain TGM
if isfield(config,'normalize')
    normalize = config.normalize;                         % Normalize empirical and shuffled TGMs by the mean and std of shuffled TGMs
else
    normalize = 'yes';
end

% Set up recurrence range over which stats will be applied
if isfield(config,'statsrange')
    statsrange = config.statsrange(1):config.statsrange(end);
else
    statsrange = 1:30;
end

% Adjust to be a factor of warped frequency in case of brain time ref dimension
if strcmp(refdimension.dim,'braintime')
    statsrange = statsrange/warpfreq;
end

%% statistically test TGM
% FIRST LEVEL PERMUTATION
% % Pre-allocate
modepow_shuff = zeros(1,numel(statsrange));
fullspec_shuff = zeros(1,numel(statsrange));
permTGM = zeros(numperms1,size(TGM,1),size(TGM,2));

% First level permutations
for perm1 = 1:numperms1
    fprintf('First level permutation number %i\n', perm1);
    clabel = clabel(randperm(numel(clabel)));
    [permTGM(perm1,:,:),~] = mv_classify_timextime(cfg_mv, bt_data.trial, clabel);
end

% If normalize, calculate mean and std and correct
if strcmp(normalize,'yes')
    mn_permTGM = mean(permTGM,1);
    mn = mean(mn_permTGM(:));
    sd = std(mn_permTGM(:));
    
    for perm1 = 1:numperms1 % Normalize shuffled data
        permTGM(perm1,:,:) = (permTGM(perm1,:,:)-mn)./sd;
    end
    
    TGM = (TGM-mn)./sd; % Normalize empirical data
end

% Analyze first level permutation
for perm1 = 1:numperms1
    
    % Calculate autocorrelation map (AC)
    ac=autocorr2d(squeeze(permTGM(perm1,:,:)));
    
    % Run FFT over all rows and columns of the AC map
    nvecs=numel(ac(:,1));
    
    % Perform FFT over one row to get f and find out statsrange indices
    if perm1 == 1
        [~,f]=Powspek(ac(1,:),nvecs/refdimension.value);
        l = nearest(f,statsrange(1)); %minimum frequency to be tested
        h = nearest(f,statsrange(end)); %maximum frequency to be tested
        srange = l:h;
    end
    
    for vec=1:nvecs
        %1st dimenssion
        [PS,f]=Powspek(ac(vec,:),nvecs/refdimension.value);
        PS1(vec,:) = PS(srange);
        
        %2nd dimension
        [PS,f]=Powspek(ac(:,vec),nvecs/refdimension.value);
        PS2(vec,:) = PS(srange);
        
    end
    avg_PS = mean(PS1,1)+mean(PS2,1); %Mean power spectra
    modepow_shuff(perm1) = avg_PS(nearest(f,modefreq)); %Mean power spectra at mode freq
    fullspec_shuff(perm1,:) = avg_PS;
end

mean_modepow_shuff = mean(modepow_shuff);
f=f(l:h); %filter frequency vector based on range of interest

% EMPIRICAL DATA
% Calculate autocorrelation map (AC)
ac=autocorr2d(TGM);

% Size of all rows and columns
nvecs=numel(ac(:,1));

% Pre-allocate
PS1 = zeros(nvecs,numel(srange));
PS2 = zeros(nvecs,numel(srange));

% Run FFT over all rows and columns of the AC map
for vec=1:nvecs
    %1st dimenssion
    [PS,f]=Powspek(ac(vec,:),nvecs/refdimension.value);
    PS1(vec,:) = PS(srange);
    
    %2nd dimension
    [PS,f]=Powspek(ac(:,vec),nvecs/refdimension.value);
    PS2(vec,:) = PS(srange);
end
f=f(l:h); %filter frequency vector based on range of interest
avg_PS = mean(PS1,1)+mean(PS2,1); %Mean power spectra
modepow_emp = avg_PS(nearest(f,modefreq)); %Mean power at mode freq
fullspec_emp = avg_PS;

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
    createCI = true;
    if numperms1 <100
        warning('on');
        warning('No confidence interval will be displayed as the number of first level permutations is too low (<100)')
        createCI = false;
    end
    
    if createCI == true
        for f_ind = 1:numel(f)
            % Confidence interval
            low_CI(f_ind,:) = prctile(fullspec_shuff(:,f_ind),2.5);
            hi_CI(f_ind,:) = prctile(fullspec_shuff(:,f_ind),97.5);
        end
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
    p1 = plot(f,fullspec_emp,'LineWidth',2,'Color','b');
    
    % Only create confidence intervals with more than 100 permutations
    if createCI == true
    plot(f,low_CI,'LineStyle','-','LineWidth',0.5,'Color','k');
    plot(f,hi_CI,'LineStyle','-','LineWidth',0.5,'Color','k');
    patch([f fliplr(f)],[low_CI' fliplr(hi_CI')], 1,'FaceColor', 'black', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    end
    
    if strcmp(refdimension.dim,'braintime') %warp freq line is dependent on clock (warped freq) or brain time (1 hz)
        p2 = line([1 1], [0 max(fullspec_emp)],'color',[1 0 1],'LineWidth',4); %Line at warped freq
        xlabel('Recurrence frequency (factor of warped freq)')
    else
        p2 = line([warpfreq warpfreq], [0 max(fullspec_emp)],'color',[1 0 1],'LineWidth',4); %Line at warped freq
        xlabel('Recurrence frequency')
    end
    p2.Color(4) = 0.45;
    
    % set up axes
    ylabel('Mean amplitude')
    title(sprintf('Amplitude across recurrence rates'))
    legend([p1 p2],{'Empirical amplitude','Warped frequency'});
end

%% Create output structure
stats1.f = f;                                         % Frequency vector
stats1.empTGM = TGM;                                  % Empirical TGM
stats1.shuffTGM = permTGM;                            % Nperm1 shuffled TGMs
stats1.empspec = fullspec_emp;                        % Power spectrum of average empirical data
stats1.shuffspec = fullspec_shuff;                    % Power spectrum of average permutation data
stats1.empmode = modepow_emp;                         % Power at mode freq of epirical data
stats1.shuffmode = modepow_shuff;                     % Power at mode freq of permutation data