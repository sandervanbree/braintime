function [stats2] = bt_statslevel2(config,config_clus,stats1)
% Performs second level statistics of periodicity power spectra, and tests
% which TGM or diagonal datapoints are significantly higher than the
% permutation distribution. Both the empirical and permuted spectra are
% z-scored.
%
% Second level statistics uses the first level permuted data.
% Specifically, nperm2 times, grab a random power spectrum of every
% participants' nperm1 pool and average it. Then, for each frequency, the
% p-value is calculated as the percentile of the empirical power in the
% nperm2 null distribution.
%
% Finally, performance in the TGM or diagonal is tested for statistical 
% significance using cluster correction. For this, the toolbox calls MVPA
% Light.
%
% Use:
% [stats2] = bt_statslevel2(config, config_clus, stats1), where stats1 has
% one cell for each participant.
%
% Input Arguments:
% config
%   - numperms2      % Number of permutations on the second statistical
%                    % level.
%                    %
%   - multiplecorr   % Correction method for multiple comparisons.
%                    % 'none' (default), 'fdr' for adjusted p-values based
%                    % on the False Discovery Rate, and 'bonferroni' for
%                    % Bonferroni correction.
%                    %
%   - cfg_clus       % Separate structure with cluster correction
%                    % parameters, used as input to MVPA Light.
%                    % For detailed information on all parameters, call
%                    % "help mv_statistics" or see ".../MVPA-Light-master/
%                    % examples/understanding_statistics.m"
%                    %
%                    % Primary parameters include:
%                    % .test              ['binomial', 'permutation']
%                    % .correctm          ['none', 'bonferroni', 'cluster']
%                    % .tail              [-1 or +1]
%                    % .n_permutations    [1000]
%                    % .clusterstatistic  ['maxsum', 'maxsize']
%                    % .alpha             [0.05]
%                    % .statistic         ['ttest' 'wilcoxon' 'mean']
%                    % .null              [0.5]
%                    % .clustercritval    [z-value; 1.96]
%                    %
% stats1             % Output structure obtained using bt_statslevel1.
%                    % Each cell contains one participant's empirical and
%                    % permuted TGMs/diagonals, the empirical and permuted
%                    % periodicity power spectra, and the tested
%                    % frequencies.
%                    %
% Output:            % 
% stats2             % Output structure that contains periodicity p-values
%                    % and tested frequencies, the employed analysis type,
%                    % and a cluster statistics structure from MVPA Light.


%% Crop away empty participant cells
stats1 = stats1(~cellfun('isempty',stats1));

%% Get information
numsubj = numel(stats1);                             % Number of participants
numperms1 = size(stats1{1}.permspec,1);              % Number of first level permutations
numperms2 = config.numperms2;                        % Number of second level permutations
mv_cfg = stats1{1}.mv_results.cfg;                   % MVPA Light results structure

refdimension = stats1{1}.refdimension;               % Find out the reference dimension (clock or brain time)

% Check
if numperms2 <10^5
    numperms2 = 10^5;
    warning(['Due to high variance in this statistical procedure, the number of second level permutations has been increased from ', num2str(numperms2),' to 100000.']);
end

%% STEP 1: Z-SCORE EMPIRICAL AND PERMUTED PERIODICITY POWER SPECTRA
for i = 1:numel(stats1)                              % For each participant...
    
    % Empirical
    mn_emp = mean(stats1{i}.empspec);
    sd_emp = std(stats1{i}.empspec);
    stats1{i}.empspec = (stats1{i}.empspec-mn_emp)./sd_emp;
    
    % Permuted
    for p = 1:size(stats1{i}.permspec,1)             % Loop through permutations
    mn_perm = mean(stats1{i}.permspec(p,:));
    sd_perm = std(stats1{i}.permspec(p,:));
    stats1{i}.permspec(p,:) = (stats1{i}.permspec(p,:)-mn_perm)./sd_perm;
    end
end

%% STEP 2: SECOND LEVEL STATISTICS OF PERIODICITY POWER SPECTRA
%% Scale the results to be of the same frequency range and length
% Determine frequency ranges of participants
for i = 1:numel(stats1)                              % Find each participant's min and max freq
    minf(i)    = min(stats1{i}.f);
    maxf(i)    = max(stats1{i}.f);
end
frange = [max(minf),min(maxf)];                      % What freq do all participants have in common?

% Filter stats1 based on common frequencies across participants
for i = 1:numel(stats1)
    st = nearest(stats1{i}.f,frange(1));             % Find nearest starting freq for each participant
    en = nearest(stats1{i}.f,frange(2));             % Find nearest ending freq for each participant
    
    ffix = @(x) x(st:en);                            % Filter all cells based on common frequencies
    stats1{i}.f         = ffix(stats1{i}.f);
    nfbins(i)           = numel(stats1{i}.f);
    stats1{i}.empspec   = ffix(stats1{i}.empspec);
    
    for i2 = 1:size(stats1{i}.permspec,1)           % Do that for each permutation
        temp(i2,:) = ffix(stats1{i}.permspec(i2,:));
    end
    stats1{i}.permspec = temp;                % And replace
    clear temp
end

% Resize stats1 cells to the same length
res = @(x) imresize(x,[1 mode(nfbins)]);             % Resize to the most common number of frequency bins to minimize data interpolation
for i = 1:numel(stats1)
    stats1{i}.f         = res(stats1{i}.f);
    fvec(i,:)           = stats1{i}.f;
    stats1{i}.empspec   = res(stats1{i}.empspec);
    
    for i2 = 1:size(stats1{i}.permspec,1)          % Do that for each permutation
        temp2(i2,:) = res(stats1{i}.permspec(i2,:));
    end
    stats1{i}.permspec = temp2;                     % And replace
end

% Create new f - averaging slight variations
f = nanmean(fvec,1);

%% Get average periodicity power spectrum across participants
% Create full empirical power spectrum
for subj = 1:numsubj
    PS_emp(subj,:) = stats1{subj}.empspec;
end
PS_emp_avg = mean(PS_emp,1);

%% Second level statistics
% Pre-allocate
perm1PS = zeros(numsubj,mode(nfbins));
perm2PS = zeros(numperms2,mode(nfbins));
progbar = round(linspace(0,numperms2,100));

% Loop through numperms2
for perm2=1:numperms2
    if ismember(perm2,progbar) % Print progress
        disp(strcat((num2str(round((perm2/numperms2)*100))),'% of second level permutations completed (for second-level stats)'));
    end
    for subj = 1:numsubj
        idx=randperm(numperms1,1); %randomly grab with replacement
        perm1PS(subj,:) = stats1{subj}.permspec(idx,:);
    end
    perm2PS(perm2,:) = mean(perm1PS,1);
end

% Mean across all 2nd level permutations (for plotting)
perms2PS_avg = mean(perm2PS,1);

%% Calculate p-values and confidence intervals
% For each frequency, find out the percentile of the empirical power in the permuted distribution
for f_ind = 1:numel(f)
    
    % p-value
    pval(f_ind) = numel(find(perm2PS(:,f_ind)>=PS_emp_avg(f_ind)))/numperms2;
    
    % Confidence interval
    low_CI(f_ind,:) = prctile(perm2PS(:,f_ind),2.5);
    hi_CI(f_ind,:) = prctile(perm2PS(:,f_ind),97.5);
end

%% Get the indices of warping frequency harmonics
if strcmp(refdimension.dim,'braintime')
    %Find the warped frequency (1 Hz)
    wfreq_i = nearest(f,1);
    
    % Find harmonics
    if f(1) - 0.5 < 0.1 % Is half the warped frequency in the tested range?
        wfreq_half_i = nearest(f,0.5); %Find half the warped frequency (0.5 Hz)
    end
    
    if f(end) - 2 > -0.1 % Is double the warped frequency in the tested range?
        wfreq_double_i = nearest(f,2); %Find double the warped frequency (2 Hz)
    end
end

%% Multiple comparisons correction
pval_corr = pval;
if isfield(config,'multiplecorr')
    if strcmp(config.multiplecorr,'fdr')
        [~,~,~,pval_corr] = fdr_bh(pval,0.05,'dep');
    elseif strcmp(config.multiplecorr,'bonferroni')
        pval_corr = pval*numel(pval);
    elseif strcmp(config.multiplecorr,'none')
        pval_corr = pval_corr;
    else
        error(['Multiple correction option "',config.multiplecorr,'" is not recognized. Please avoid capital letters or see help bt_statslevel2 for all options.'])
    end
end

% For brain time referenced data, exempt the warped frequency and its harmonics from multiple testing correction
if strcmp(refdimension.dim,'braintime')
    
    % Exempt harmonics and add harmonics in output structure
    if exist('wfreq_half_i','var') == 1 &&  exist('wfreq_double_i','var') == 1 %Exempt tested harmonics from p-value correction
        pval_corr([wfreq_half_i,wfreq_i,wfreq_double_i]) = pval([wfreq_half_i,wfreq_i,wfreq_double_i]);
        stats2.warpedfreqs = f([wfreq_half_i,wfreq_i,wfreq_double_i]);
    elseif exist('wfreq_half_i','var') == 1 &&  exist('wfreq_double_i','var') == 0
        pval_corr([wfreq_half_i,wfreq_i]) = pval([wfreq_half_i,wfreq_i]);
        stats2.warpedfreqs = f([wfreq_half_i,wfreq_i]);
    elseif exist('wfreq_half_i','var') == 0 &&  exist('wfreq_double_i','var') == 1
        pval_corr([wfreq_i,wfreq_double_i]) = pval([wfreq_i,wfreq_double_i]);
        stats2.warpedfreqs = f([wfreq_i,wfreq_double_i]);
    elseif exist('wfreq_half_i','var') == 0 &&  exist('wfreq_double_i','var') == 0
        pval_corr(wfreq_i) = pval(wfreq_i);
        stats2.warpedfreqs = f(wfreq_i);
    end
    
    fprintf(1, '\n');
    disp('##############################################################################%##################################')
    disp('The p-value at the warped frequency and harmonics (0.5, 1, and 2Hz) are exempted from multiple testing correction');
    disp('#################################################################################################################')
    fprintf(1, '\n');
end

pval_corr(pval_corr>1)=1; % Prevent >1 p values.

% Get the negative logarithm for visualization purposes
logpval = -log10(pval_corr);
if any(isinf(logpval))
    fprintf(1, '\n');
    disp('The empirical power of at least one frequency is higher than all permuted power values (yielding p = 0).');
    disp('Consider increasing the number of 2nd level permutations');
    fprintf(1, '\n');
    logpval(isinf(logpval)) = 4; %cap logpval on 4.
end

%% Plot results
% Relationship empirical and empirical power at all freqs
figure;hold on;
bt_figure('clocktime_per');

if strcmp(refdimension.dim,'braintime') % Only for brain time, separate warped frequency
    subplot(1,10,1:3);
    bt_figure('braintime_per');
    
    % Plot the warped freq
    wfreq_emp = PS_emp_avg(:,wfreq_i); % Empirical data at warped freq
    wfreq_perm = perm2PS(:,wfreq_i); % Permuted data at warped freq
    plot(2,wfreq_emp,'o','MarkerSize',9,'MarkerEdgeColor',bt_colorscheme('per_ps_emp'),'MarkerFaceColor',bt_colorscheme('per_ps_emp'));hold on; %Plot marker of empirical power
    
    if exist('wfreq_half_i','var') == 1 %Plot half the warped freq?
        half_wfreq_emp = PS_emp_avg(:,wfreq_half_i);
        half_wfreq_perm = perm2PS(:,wfreq_half_i);
        plot(1,half_wfreq_emp,'o','MarkerSize',9,'MarkerEdgeColor',bt_colorscheme('per_ps_emp'),'MarkerFaceColor',bt_colorscheme('per_ps_emp'));hold on; %Plot marker of empirical power
    else
        half_wfreq_emp = zeros(size(wfreq_emp));
        half_wfreq_perm = zeros(size(wfreq_perm));
    end
    
    if exist('wfreq_double_i','var') == 1 %Plot double the warped freq?
        double_wfreq_emp = PS_emp_avg(:,wfreq_double_i);
        double_wfreq_perm = perm2PS(:,wfreq_double_i);
        plot(3,double_wfreq_emp,'o','MarkerSize',9,'MarkerEdgeColor',bt_colorscheme('per_ps_emp'),'MarkerFaceColor',bt_colorscheme('per_ps_emp'));hold on; %Plot marker of empirical power
    else
        double_wfreq_emp = zeros(size(wfreq_emp));
        double_wfreq_perm = zeros(size(wfreq_perm));
    end
    
    % Put all the three frequencies in one matrix
    allPerm = [half_wfreq_perm,wfreq_perm,double_wfreq_perm];
    allEmp = [half_wfreq_emp,wfreq_emp,double_wfreq_emp];
    
    % Create a Violin plot. If this does not work because of an older MATLAB version, make a boxplot instead
    try
        violinplot(allPerm,'test','ShowData',false,'ViolinColor',bt_colorscheme('confidenceinterval'),'MedianColor',bt_colorscheme('per_ps_perm'),'BoxColor',[0.7 0.7 0.7],'EdgeColor',[1 1 1],'ViolinAlpha',0.15);
        % Set legend
        h = get(gca,'Children');
        l2 = legend(h([numel(h) 3]),'Empirical power','Permuted power');
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
        l2 = legend([h(2),toMark(1)],'Empirical power','Permuted power');
        set(l2,'Location','best');
    end
    
    % Set up axes
    ylabel('Average power (Z)');
    title('2nd level periodicity')
    
    % Adapt xticklabels depending on which freqs are in range
    if sum(half_wfreq_emp) == 0 && sum(double_wfreq_emp) == 0
        xticklabels({'0.5 Hz (UNTESTED)','1 Hz (warp)','2 Hz (UNTESTED)'});
    elseif sum(half_wfreq_emp) ~= 0 && sum(double_wfreq_emp) == 0
        xticklabels({'0.5 Hz','1 Hz (warp)','2 Hz (UNTESTED)'});
    elseif sum(half_wfreq_emp) == 0 && sum(double_wfreq_emp) ~= 0
        xticklabels({'0.5 Hz (UNTESTED)','1 Hz (warp)','2 Hz'});
    elseif sum(half_wfreq_emp) ~= 0 && sum(double_wfreq_emp) ~= 0
        xticklabels({'0.5 Hz','1 Hz (warp)','2 Hz'});
    end
    
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    l2.FontSize = bt_plotparams('FontSizeLegend');
    
    title('2nd level periodicity')
    
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    
    % Prepare new subplot
    subplot(1,10,5:10)
    hold on
end

yyaxis left
p1 = plot(f,PS_emp_avg,'LineStyle','-','LineWidth',3,'Color',bt_colorscheme('per_ps_emp'),'Marker','none'); %Mean across 2nd level perms
p2 = plot(f,perms2PS_avg,'LineStyle','-','LineWidth',2,'Color',bt_colorscheme('per_ps_perm'),'Marker','none'); %Mean across 2nd level perms
xlabel('Periodicity frequency (Hz)')
ylabel('Average power (Z)')

% Plot confidence interval
c2 = plot(f,low_CI,'LineWidth',0.5,'Color',[0 0 0],'Marker','none','LineStyle','none');
c3 = plot(f,hi_CI,'LineWidth',0.5,'Color',[0 0 0],'Marker','none','LineStyle','none');
p3 = patch([f fliplr(f)],[low_CI' fliplr(hi_CI')], 1,'FaceColor', bt_colorscheme('confidenceinterval'), 'EdgeColor', 'none', 'FaceAlpha', 0.15);

% p-value axis
yyaxis right
ax = gca;
ax.YAxis(1).Color = 'black';
p4 = plot(f,logpval,'LineStyle','-','LineWidth',2,'Color',bt_colorscheme('pval'));
p4.Color(4) = 0.3;
ylim([-0.05,4])
ylabel('-log10 p-value')

% plot star at every significant frequency
yyaxis left
ax = gca;
ax.YAxis(2).Color = bt_colorscheme('pval');
sigind = find(pval_corr<=0.05);
if isempty(sigind) ~= 1
    p5 = plot(f(sigind),0,'r*','MarkerSize',10,'LineWidth',1.5,'color',bt_colorscheme('sigstar'));
    % legend
    l2 = legend([p1 p2 p3 p4 p5(1)],{'Empirical','Permuted', 'Conf. interv.', '-log10 p-value', 'p <= 0.05'});
else
    l2 = legend([p1 p2 p3 p4],{'Empirical','Permuted', 'Conf. interv.', '-log10 p-value'});
end
set(l2,'Location','best')


% add title
title('Periodicity power spectra (2nd level stats)');

% Adapt font
set(gca,'FontName',bt_plotparams('FontName'));
set(gca,'FontSize',bt_plotparams('FontSize'));
l2.FontSize = bt_plotparams('FontSizeLegend');

%% Adapt output variable
% Add freq and pval information to output
stats2.periodicity.pval_corr = pval_corr;
stats2.periodicity.frequency = f;
stats2.analysistype = ['Periodicity was measured across the ', upper(stats1{1}.method)];

% Print results for brain time
disp('See the output structure for:')
disp('(1) the (corrected) p-values per frequency');
disp('(2) the range of tested frequencies');
disp('(3) for brain time, the (corrected) p-value of 0.5x, 1x, and 2x the warped frequency');
disp('(4) The analysis type employed');

%% STEP 3: TEST SIGNIFICANCE OF TGM OR DIAGONAL USING MVPA LIGHT

for i = 1:numsubj                                              % Prepare group structure
    mv_clus{i} = stats1{i}.mv_results;
end

config_clus.metric = mv_cfg.metric;                            % Append config parameters that are fixed
config_clus.design = 'within';

try clus_level2 = mv_statistics(config_clus, mv_clus);         % Run cluster statistics
catch
    clus_level2.p = ['No clusters reached statistical significance '...
        'cfg_clus.clustercritval may be set too conservatively.'];
    clus_level2.mask = [];
end

mapmask = clus_level2.mask;                                    % Extract mask
mv_clus_average = mv_combine_results(mv_clus, 'average');      % Combine MVPA results

% Plotting
if strcmp(mv_clus_average.function,'mv_classify_timextime')    % For cross temporal generalization
    
    TGM_avg = mv_clus_average.perf{1};                         % Fetch TGM average
    
    figure;hold on;bt_figure('cluster');
    cfg_plot = [];
    cfg_plot.ylim = [1 size(TGM_avg,1)];
    cfg_plot.xlim = [1 size(TGM_avg,2)];
    cfg_plot.colorbar = 1;
    cfg_plot.colorbar_location = 'EastOutside';
    cfg_plot.colorbar_title = mv_cfg.metric;
    cfg_plot.climzero = config_clus.null;
    mv_plot_2D(cfg_plot,TGM_avg);hold on;
    axis square;
    ax = gca;
    ax.YLim = [1 size(TGM_avg,1)];
    ax.XLim = [1 size(TGM_avg,2)];
    colormap(bt_colorscheme('TGM'));freezeColors;
    xlabel('Test data (bin)')
    ylabel('Training data (bin)')
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    title('Average TGM with significant clusters');
    
    if isempty(mapmask)~=1
        % Add contour
        contour(1:size(TGM_avg,1),1:size(TGM_avg,2),mapmask,1,'linewidth', 2, 'color',bt_colorscheme('sigcluster_TGM'));hold on;
    else
        warning('No significant clusters were found. Try changing cfg_clust.clustercritval.')
    end
    
elseif strcmp(mv_clus_average.function,'mv_classify_across_time')    % For diagonal classification
    
    
    time_plot = 1:size(stats1{1}.mv_results.perf,1);     % in data bins because time vectors may differ for brain time
    result_average = mv_combine_results(mv_clus, 'average');
    
    figure;hold on;bt_figure('cluster');                 % Plot average results
    mv_plot_result(result_average, time_plot, 'mask', clus_level2.mask, 'new_figure', 0)
    h = findobj(gca,'Type','line');
    for i = numel(h)
        h(i).LineWidth = h(i).LineWidth+1;
    end
    
    xlabel('Test data (bin)')
    ylabel(mv_cfg.metric);
    title('Average empirical Diagonal (classifier time course)');
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    
end

%% Adapt output variable
% Add cluster correction information to output (from MVPA Light)
stats2.clusters = clus_level2;
end