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
%                    % TGM. Critical: a difference in permuted and
%                    % empirical analyses should not be caused by
%                    % differences in classification parameters.
%                    %
%   - numperms1      % Number of permutations on the first statistical
%                    % level.
%                    %
%   - recurrencefoi  % Range of recurrence rates to be statistically tested
%                    % in the TGM.
%                    %
%   - figure         % 'yes' (default) or 'no': display statistical results
%                    %
% bt_data            % Data structure obtained from bt_TGMquant. Contains:
%                    % TGM or its AC map, quantification of it, and config
%                    % details saved for later retrieval.
%                    %
% bt_TGMquant        % Recurrence power spectrum, TGM, and mapmethod used
%                    % to quantify recurrence.
%                    %
% Output:            %
% stats1             % Output structure which contains power spectra of
%                    % the empirical data, permutation data, and the
%                    % associated frequency vector.

%% Get information
numperms1 = config.numperms1;                             % Number of first level permutations
warpfreq = bt_TGMquant.warpfreq;                          % Warped frequency (frequency of the carrier)
TGM = bt_TGMquant.TGM;                                    % Time Generalization Matrix of the data
refdimension = bt_TGMquant.refdimension;                  % Reference dimension used
recurrencefoi = bt_TGMquant.recurrencefoi;                % Range of tested TGM recurrence frequencies
mapmethod = bt_TGMquant.mapmethod;                        % Check whether analysis is done over TGM or AC map
pspec_emp = bt_TGMquant.pspec_emp;                        % Recurrence power spectrum of empirical data
clabel = config.clabel;                                   % Classification labels
cfg_mv = config.mvpacfg;                                  % MVPA Light configuration structured used to obtain TGM

if isfield(config,'normalize')
    normalize = config.normalize;                         % Normalize empirical and permuted TGMs by the mean and std of permuted TGMs
else
    normalize = 'yes';
end

% Set up recurrence range over which stats will be applied
powspecrange = recurrencefoi;

%% statistically test TGM
% FIRST LEVEL PERMUTATION
% % Pre-allocate
pspec_perm = zeros(1,numel(powspecrange));
permTGM = zeros(numperms1,size(TGM,1),size(TGM,2));

% First level permutations
for perm1 = 1:numperms1
    fprintf('First level permutation number %i\n', perm1);
    clabel = clabel(randperm(numel(clabel)));
    [permTGM(perm1,:,:),~] = mv_classify_timextime(cfg_mv, bt_data.trial, clabel);
end

% Analyze first level permutation
for perm1 = 1:numperms1
    
    % Perform analysis over TGM or its autocorrelation map (AC)?
    if strcmp(mapmethod,'tgm')
        mp = squeeze(permTGM(perm1,:,:));
    else
        mp = autocorr2d(squeeze(permTGM(perm1,:,:)));
    end
    
    % Run FFT over all rows and columns of the AC map
    nrows = numel(mp(:,1));
    ncols = numel(mp(:,2));
        
    for row = 1:nrows % Perform FFT over rows
        
        if row == 1 % For the first row, perform a test analysis
            [~,f]=Powspek(mp(1,:),nrows/refdimension.value);
            l = nearest(f,powspecrange(1)); %minimum frequency to be tested
            h = nearest(f,powspecrange(end)); %maximum frequency to be tested
            ps_range = l:h; % this is the range of frequencies desired
        end
        
        % 1st dimension
        [PS,f]=Powspek(mp(row,:),nrows/refdimension.value);
        PS1(row,:) = PS(ps_range); % restrict do desired range
    end
    
    for col = 1:ncols % Perform FFT over columns
        
        if col == 1 % For the first column, perform a test analysis
            [~,f]=Powspek(mp(1,:),ncols/refdimension.value);
            l = nearest(f,powspecrange(1)); %minimum frequency to be tested
            h = nearest(f,powspecrange(end)); %maximum frequency to be tested
            ps_range = l:h; % this is the range of frequencies desired
        end
        
        % 2nd dimension
        [PS,f]=Powspek(mp(:,col),ncols/refdimension.value);
        PS2(col,:) = PS(ps_range); % restrict do desired range
    end
    
    avg_PS = mean(PS1,1)+mean(PS2,1); %Mean power spectra
    pspec_perm(perm1,:) = avg_PS;
end
    
f=f(ps_range); %filter frequency vector based on range of interest

% Only calculate confidence interval and plot stats if desired
if isfield(config,'figure')
    if strcmp(config.figure,'no')
        figopt = 0;
    else
        figopt = 1; %Default yes   
    end
end

if figopt == 1
    %% Create confidence interval for each frequency bin
    createCI = true;
    if numperms1 <20
        warning('on');
        warning('No confidence interval will be displayed as the number of first level permutations is too low (<20)')
        createCI = false;
    end
    
    if createCI == true
        for f_ind = 1:numel(f)
            % Confidence interval
            low_CI(f_ind,:) = prctile(pspec_perm(:,f_ind),2.5);
            hi_CI(f_ind,:) = prctile(pspec_perm(:,f_ind),97.5);
        end
    end
    
    %% Plot results
    figure;hold on;set(gcf, 'WindowState', 'maximized'); % create full screen figure

    if strcmp(refdimension.dim,'braintime') % Only for brain time, separate warped frequency
        subplot(1,10,1:3)
        wfreq = nearest(f,1); %Find the warped frequency (1 Hz)
        plot(pspec_emp(wfreq),'o','MarkerSize',6,'MarkerEdgeColor','blue','MarkerFaceColor','blue');hold on; %Plot marker of empirical power
               
        % Create a Violin plot. If this does not work because of an older MATLAB version, make a boxplot instead
        try
            violinplot(pspec_perm(:,wfreq),'test','ShowData',false,'ViolinColor',[0.8 0.8 0.8],'MedianColor',[0 0 0],'BoxColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'ViolinAlpha',0.8);
            % Set legend
            h = get(gca,'Children');
            l2 = legend(h([9 3]),'Empirical (emp) recurrence power','Permuted (perm) recurrence power');
            set(l2,'Location','best');
        catch
            boxplot(pspec_perm(:,wfreq))           
            % Set up y-axis
            maxy = max([pspec_perm(:,wfreq)',pspec_emp(wfreq)]);
            ylim([min(pspec_perm(:,wfreq))*0.9,maxy*1.1]); % slightly below min and max            
            % Set legend
            h = get(gca,'Children');
            h(1).Children(1).Color = [1 1 1];
            h(1).Children(2).LineWidth = 3;
            h(1).Children(2).Color = [0 0 0];
            h(1).Children(2).LineWidth = 6;
            h(1).Children(3).Color = [0.6 0.6 0.6];
            h(1).Children(3).LineWidth = 3;
            h(1).Children(4).Color = [0.6 0.6 0.6];
            h(1).Children(4).LineWidth = 3;
            h(1).Children(5).Color = [0.6 0.6 0.6];
            h(1).Children(5).LineWidth = 3;
            h(1).Children(6).Color = [0.6 0.6 0.6];
            h(1).Children(6).LineWidth = 3;
            h(1).Children(7).Color = [0.6 0.6 0.6];
            
            l2 = legend([h(2),h(1).Children(2)],'Empirical (emp) recurrence power','Permuted (perm) recurrence power');
            set(l2,'Location','best');
        end
        
        % Set up axes
        ylabel('Recurrence power');
        xticklabels(' ');
        xlabel('Warped frequency (1 Hz)');
        title('1st level recurrence at warped frequency (1Hz)')
        
        % Adapt font
        set(gca,'FontName','Arial')
        set(gca,'FontSize',16)
        
        % Now plot recurrence power spectrum
        subplot(1,10,5:10)
        hold on
    end
    
    yyaxis left
    p1 = plot(f,pspec_emp,'LineStyle','-','LineWidth',3,'Color','b'); %Mean across 1st level perms
    p2 = plot(f,mean(pspec_perm,1),'LineStyle','-','LineWidth',2,'Color',[0.3 0.3 0.3]); %Mean across 1st level perms
    xlabel('Recurrence frequency')
    ylabel('Mean power across participants')
    
    if strcmp(refdimension.dim,'braintime') %warp freq line is dependent on clock (warped freq) or brain time (1 hz)
        p3 = line([1 1], [0 max(pspec_emp)],'color',[1 0 1],'LineWidth',4); %Line at warped freq
        xlabel('Recurrence frequency (factor of warped freq)')
    else
        p3 = line([warpfreq warpfreq], [0 max(pspec_emp)],'color',[1 0 1],'LineWidth',4); %Line at warped freq
        xlabel('Recurrence frequency')
    end
    p3.Color(4) = 0.45;
    
    % Plot confidence interval
    if createCI == true
        c2 = plot(f,low_CI,'LineStyle','-','LineWidth',0.5,'Color','k');
        c3 = plot(f,hi_CI,'LineStyle','-','LineWidth',0.5,'Color','k');
        p4 = patch([f fliplr(f)],[low_CI' fliplr(hi_CI')], 1,'FaceColor', 'black', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        
        % legend
        l2 = legend([p1 p2 p3 p4],{'Average emp spectrum','Average perm spectrum', 'Warped frequency', 'Conf. interv. perm spectrum'});
    else
        l2 = legend([p1 p2 p3],{'Average emp spectrum','Average perm spectrum','Warped frequency'});
    end
    set(l2,'Location','best')
    
    % add title
    title('Recurrence power spectra (1st level stats)');
    
    % Adapt font
    set(gca,'FontName','Arial')
    set(gca,'FontSize',16)
    
    % Notify user about lack of p-values
    disp('p-values are calculated in the second-level statistics (bt_TGMstatslevel2)');
end

%% Create output structure
stats1.f = f;                                         % Frequency vector
stats1.empTGM = TGM;                                  % Empirical TGM
stats1.permTGM = permTGM;                             % Nperm1 permuted TGMs
stats1.empspec = pspec_emp;                           % Power spectrum of average empirical data
stats1.permspec = pspec_perm;                         % Power spectrum of average permutation data
stats1.refdimension = refdimension;                   % Save reference dimension (clock or brain time)