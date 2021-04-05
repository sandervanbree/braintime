function [stats1] = bt_statslevel1(config, bt_data, quant)
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
%                    %
%   - numperms1      % Number of permutations on the first statistical
%                    % level.
%                    %
%   - periodicityfoi % Range of periodicity rates to be statistically tested
%                    % in the TGM.
%                    %
%   - figure         % 'yes' (default) or 'no': display statistical results
%                    %
% bt_data            % Data structure obtained from bt_TGMquant. Contains:
%                    % TGM or its AC map, quantification of it, and config
%                    % details saved for later retrieval.
%                    %
% bt_TGMquant        % MVPA configuration structure, periodicity power
%                    % spectrum, TGM, and mapmethod used to quantify
%                    % periodicity.
%                    %
% Output:            %
% stats1             % Output structure which contains power spectra of
%                    % the empirical data, permutation data, and the
%                    % associated frequency vector.

%% Get information
numperms1 = config.numperms1;                             % Number of first level permutations
MVPAcfg = quant.MVPAcfg;                            % MVPA Light configuration structured used to obtain TGM
warpfreq = quant.warpfreq;                          % Warped frequency (frequency of the carrier)
TGM = quant.TGM;                                    % Time Generalization Matrix of the data
timevec = quant.timevec;                            % Time vector (different for brain and clock time referencing)
refdimension = quant.refdimension;                  % Reference dimension used
periodicityfoi = quant.periodicityfoi;              % Range of tested TGM periodicity frequencies
mapmethod = quant.mapmethod;                        % Check whether analysis is done over TGM, AC map, or diag
pspec_emp = quant.pspec_emp;                        % Periodicity power spectrum of empirical data
clabel = config.clabel;                                   % Classification labels

% Set up periodicity range over which stats will be applied
powspecrange = periodicityfoi;

%% statistically test TGM
% FIRST LEVEL PERMUTATION
% % Pre-allocate
permTGM = zeros(numperms1,size(TGM,1),size(TGM,2));

% First level permutations
for perm1 = 1:numperms1
    fprintf('First level permutation number %i\n', perm1);
    clabel = clabel(randperm(numel(clabel)));
    [permTGM(perm1,:,:),~] = mv_classify_timextime(MVPAcfg, bt_data.trial, clabel);
end

% Analyze first level permutation
for perm1 = 1:numperms1
    
    if strcmp(mapmethod,'tgm')
        mp = squeeze(permTGM(perm1,:,:));
    elseif strcmp(mapmethod,'ac')
        mp = autocorr2d(squeeze(permTGM(perm1,:,:)));
    elseif strcmp(mapmethod,'diag')
        mp = diag(squeeze(permTGM(perm1,:,:)))';
    end

[PS,f] = bt_fft(mp,powspecrange,timevec);
pspec_perm(perm1,:) = PS;

end
    
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
    
    %% Get the indices of warping frequency harmonics 
    if strcmp(refdimension.dim,'braintime')
        %Find the warped frequency (1 Hz)
        wfreq_i = nearest(f,1);
        
        % Find harmonics
        if f(1) - 0.5 < 0.1 % Is half the warped frequency in the tested range?
        wfreq_half_i = nearest(f,0.5); %Find half the warped frequency (0.5 Hz)
        else
            warning("Periodicity may appear at 0.5x the warped frequency, but this is outside the tested range.");
        end
        
        if f(end) - 2 > -0.1 % Is double the warped frequency in the tested range?
        wfreq_double_i = nearest(f,2); %Find double the warped frequency (2 Hz)
        else
            warning("Periodicity may appear at 2x the warped frequency, but this is outside the tested range.");
        end
    else
        wfreq_i = nearest(f,warpfreq);
    end
        
    %% Plot results
    figure;hold on;bt_figure('clocktime_per');

    if strcmp(refdimension.dim,'braintime') % Only for brain time, separate warped frequency
        bt_figure('braintime_per');
        subplot(1,10,1:3)
        
        % Plot the warped freq
        wfreq_emp = pspec_emp(:,wfreq_i); % Empirical data at warped freq
        wfreq_perm = pspec_perm(:,wfreq_i); % Permuted data at warped freq
        plot(2,wfreq_emp,'o','MarkerSize',6,'MarkerEdgeColor',bt_colorscheme('per_ps_emp'),'MarkerFaceColor',bt_colorscheme('per_ps_emp'));hold on; %Plot marker of empirical power
        
        if exist('wfreq_half_i','var') == 1 %Plot half the warped freq?
        half_wfreq_emp = pspec_emp(:,wfreq_half_i);
        half_wfreq_perm = pspec_perm(:,wfreq_half_i);
        plot(1,half_wfreq_emp,'o','MarkerSize',6,'MarkerEdgeColor',bt_colorscheme('per_ps_emp'),'MarkerFaceColor',bt_colorscheme('per_ps_emp'));hold on; %Plot marker of empirical power
        else
        half_wfreq_emp = zeros(size(wfreq_emp));
        half_wfreq_perm = zeros(size(wfreq_perm));
        end
        
        if exist('wfreq_double_i','var') == 1 %Plot double the warped freq?
        double_wfreq_emp = pspec_emp(:,wfreq_double_i);
        double_wfreq_perm = pspec_perm(:,wfreq_double_i);
        plot(3,double_wfreq_emp,'o','MarkerSize',6,'MarkerEdgeColor',bt_colorscheme('per_ps_emp'),'MarkerFaceColor',bt_colorscheme('per_ps_emp'));hold on; %Plot marker of empirical power
        else 
        double_wfreq_emp = zeros(size(wfreq_emp));
        double_wfreq_perm = zeros(size(wfreq_perm));
        end
        
        % Put all the three frequencies in one matrix
        allPerm = [half_wfreq_perm,wfreq_perm,double_wfreq_perm];
        allEmp = [half_wfreq_emp,wfreq_emp,double_wfreq_emp];
        
        % Create a Violin plot. If this does not work because of an older MATLAB version, make a boxplot instead
        try
           violinplot(allPerm,'test','ShowData',false,'ViolinColor',[0.8 0.8 0.8],'MedianColor',[0 0 0],'BoxColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'ViolinAlpha',0.3);
           % Set legend
            h = get(gca,'Children');
            l2 = legend(h([numel(h) 3]),'Empirical periodicity power','Permuted periodicity power');
            set(l2,'Location','best');
        catch
            boxplot(allPerm)           
            % Set up y-axis
            maxy = max([allPerm(:);allEmp(:)]);
            miny = min([allPerm(:);allEmp(:)]);
            ylim([miny*0.9,maxy*1.1]); % slightly below min and max            
            % Set legend
            toMark = findobj('Color','red','LineStyle','-');
            h = get(gca,'Children');
            l2 = legend([h(2),toMark(1)],'Empirical periodicity power','Permuted periodicity power');
            set(l2,'Location','best');
        end
        
        % Set up labels
        ylabel('Periodicity power');
        title('1st level periodicity')
        
        % Adapt xticklabels depending on which freqs are in range
        if sum(half_wfreq_emp) == 0 && sum(double_wfreq_emp) == 0
            xticklabels({'0.5Hz (UNTESTED)','1Hz (warp)','2Hz (UNTESTED)'});
        elseif sum(half_wfreq_emp) ~= 0 && sum(double_wfreq_emp) == 0
            xticklabels({'0.5Hz','1Hz (warp)','2Hz (UNTESTED)'});
        elseif sum(half_wfreq_emp) == 0 && sum(double_wfreq_emp) ~= 0
            xticklabels({'0.5Hz (UNTESTED)','1Hz (warp)','2Hz'});
        elseif sum(half_wfreq_emp) ~= 0 && sum(double_wfreq_emp) ~= 0
            xticklabels({'0.5Hz','1Hz (warp)','2Hz'});
        end
        
        % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    l2.FontSize = bt_plotparams('FontSizeLegend');
    
        % Now plot periodicity power spectrum
        subplot(1,10,5:10)
        hold on
    end
    
    p1 = plot(f,pspec_emp,'LineStyle','-','LineWidth',3,'Color',bt_colorscheme('per_ps_emp')); %Mean across 1st level perms
    p2 = plot(f,mean(pspec_perm,1),'LineStyle','-','LineWidth',2,'Color',bt_colorscheme('per_ps_perm')); %Mean across 1st level perms
    xlim([f(1) f(end)]);
    xlabel('Periodicity frequency')
    ylabel('Power')
    
    p3 = line([f(wfreq_i) f(wfreq_i)], [0 max(pspec_emp)],'color',bt_colorscheme('warpingsignal'),'LineWidth',4); %Line at warped freq
    p3.Color(4) = 0.85;
    if strcmp(refdimension.dim,'braintime') %warp freq line is dependent on clock (warped freq) or brain time (1 hz)
        xlabel('Periodicity frequency (factor of warped freq)')
    else
        xlabel('Periodicity frequency')
    end
        
    % Plot confidence interval
    if createCI == true
        c2 = plot(f,low_CI,'LineStyle','-','LineWidth',0.5,'Color','k','LineStyle','none');
        c3 = plot(f,hi_CI,'LineStyle','-','LineWidth',0.5,'Color','k','LineStyle','none');
        p4 = patch([f fliplr(f)],[low_CI' fliplr(hi_CI')], 1,'FaceColor', bt_colorscheme('confidenceinterval'), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        
        % legend
        l2 = legend([p1 p2 p3 p4],{'Empirical spectrum','Permuted spectrum', 'Warped frequency', 'Conf. interv. perm spectrum'});
    else
        l2 = legend([p1 p2 p3],{'Empirical spectrum','Permuted spectrum','Warped frequency'});
    end
    set(l2,'Location','best')
    
    % add title
    title('Periodicity power spectra (1st level stats)');
    
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    l2.FontSize = bt_plotparams('FontSizeLegend');

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
stats1.mapmethod = mapmethod;                         % Save whether analysis is done over TGM, AC map, or diag
stats1.MVPAcfg = MVPAcfg;                             % MVPA Light configuration structured used to obtain TGM

