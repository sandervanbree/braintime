function [stats2] = bt_statslevel2(config, stats1)
% Performs second level statistics of periodicity power spectra, and tests
% which TGM datapoints are significantly higher than the distribution of
% permuted TGMs.
%
% - Second level statistics uses the first level permutated data. Numperm2
% times, grab a random power spectrum of every participants' numperm1 pool
% of permutated data and average. Then, for every frequency, the pvalue is
% the percentile of the average empirical data's periodicity power spectra
% in the null distribution of the permuted data's periodicity power spectra.
%
% - TGM significance testing compares each empirical TGM datapoint to
% the same datapoint in the permuted distribution and tests whether
% classifier performance exceeds 5%. Empirical and permuted  TGMs are
% normalized using the mean and standard deviation of the permuted TGMs.
%
% Use:
% [pval] = bt_statslevel2(config, stats1), where stats1 has one cell
% for each participant.
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
%   - cluster_p      % The threshold for significant clusters in the
%                    % cluster correction (default = 0.05).
%                    %
% stats1             % Output structure obtained using bt_statslevel1.
%                    % Each cell contains one participants' empirical TGM,
%                    % permuted TGMs, periodicity power spectrum of the
%                    % empirical and permuted TGM, and the associated
%                    % frequency vector.
%                    %
% Output:            %
% UPDATE

%% Crop away empty participant cells
stats1 = stats1(~cellfun('isempty',stats1));

%% Get information
numsubj = numel(stats1);                             % Number of participants
numperms1 = size(stats1{1}.permspec,1);              % Number of first level permutations
mapmethod = stats1{1}.mapmethod;                     % Save whether analysis is done over TGM, AC map, or diag
MVPAcfg = stats1{1}.MVPAcfg;                         % MVPA Light configuration structured used to obtain TGM

cluster_p = bt_defaultval(config,'cluster_p',0.05);          % Cluster-correction p-value threshold
cluster_n = bt_defaultval(config,'cluster_n',10);            % Cluster-correction number of clusters
cluster_smooth = bt_defaultval(config,'cluster_smooth',0);   % Smoothing for better Cluster-correction
numperms2 = bt_defaultval(config,'numperms2',10^5);          % Number of second level permutations

refdimension = stats1{1}.refdimension;               % Find out the reference dimension (clock or brain time)

% Sanity check
if numperms2 <10^5
    numperms2 = 10^5;
    warning(['Due to high variance in this statistical procedure, the number of second level permutations has been increased from ', num2str(numperms2),' to 100000.']);
end

%% STEP 1: SECOND LEVEL STATISTICS OF PERIODICITY POWER SPECTRA
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
        error(['Multiple correction option "',config.multiplecorr,'" is not recognized. Please avoid capital letters or see help bt_TGMstatslevel2 for all options.'])
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
    
    % Set up axes
    ylabel('Periodicity power');
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
xlabel('Periodicity frequency')
ylabel('Mean power across participants')

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
    l2 = legend([p1 p2 p3 p4 p5(1)],{'Empirical spectrum','Permuted spectrum', 'Conf. interv.', '-log10 p-value', 'p <= 0.05'});
else
    l2 = legend([p1 p2 p3 p4],{'Empirical spectrum','Permuted spectrum', 'Conf. interv.', '-log10 p-value'});
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
stats2.mapmethod = mapmethod;

% Print results for brain time
disp('See the output structure for:')
disp('(1) the (corrected) p-values per frequency');
disp('(2) the range of tested frequencies');
disp('(3) for brain time, the (corrected) p-value of 0.5x, 1x, and 2x the warped frequency');

%% STEP 2: TEST SIGNIFICANCE OF TGM OR DIAGONAL
% Pre-allocate
empset_TGM_Z = zeros(size(stats1{1}.empTGM,1),size(stats1{1}.empTGM,2),numel(stats1));
permset_TGM_Z = zeros(size(stats1{1}.permTGM,1),size(stats1{1}.permTGM,2),size(stats1{1}.permTGM,3),numel(stats1));
empset_TGM_ori = zeros(size(stats1{1}.empTGM,1),size(stats1{1}.empTGM,2),numel(stats1));

% Put each participant's empirical and permuted data into one matrix, whilst z-scoring based on permuted data
disp('To test for significant points in the TGM, the toolbox normalizes to');
disp('the mean and standard deviation of permuted TGMs');
for i = 1:numel(stats1)
    % Calculate mean and std of the permuted TGMs
    mn = mean(stats1{i}.permTGM(:));
    sd = std(stats1{i}.permTGM(:));
    
    % Save original
    empset_TGM_ori(:,:,i) = stats1{i}.empTGM;
    
    % Subtract the mean and divide by std to normalize across participants
    empset_TGM_Z(:,:,i) = (stats1{i}.empTGM-mn)./sd;
    permset_TGM_Z(:,:,:,i) = (stats1{i}.permTGM-mn)./sd;
end

% For every participant, average the permuted TGMs
permset_TGM_Z = squeeze(mean(permset_TGM_Z,1));

% Compute average TGM
TGMavg = squeeze(mean(empset_TGM_ori,3));

% Determine number of cluster permutations based on TGM size (larger TGM =
% fewer permutations).
% TGMsz = numel(stats1{1}.empTGM(:));
% cluster_nperm = 1./(TGMsz/(10^6));
% cluster_nperm = round(cluster_nperm.*10^4);
% if cluster_nperm < 10^4 % Minimum
%     cluster_nperm = 10^4;
% elseif cluster_nperm > 10^5 % Maximum
%     cluster_nperm = 10^5;
% end

cluster_nperm = 10^4;

% Filter diagonal if diagonal is analyzed
if strcmp(mapmethod,'diag')
    for i = 1:size(empset_TGM_Z,3)
        emp_clus(:,i) = diag(empset_TGM_Z(:,:,i))';
        perm_clus(:,i) = diag(permset_TGM_Z(:,:,i))';
    end
else % Else just keep the whole TGM
    emp_clus = empset_TGM_Z;
    perm_clus = permset_TGM_Z;
end

% Smoothing
if cluster_smooth~=0
    emp_clus = imgaussfilt(emp_clus,cluster_smooth);
    perm_clus = imgaussfilt(perm_clus,cluster_smooth);
end

% Notify
disp('Testing for significant clusters, this may take a while...');

% Perform cluster correction (Maris & Oostenveld, 2007; J Neurosci Methods)
% where the empirical TGMs are compared against the average permuted TGMs
[c_TGM,p_TGM,~,~] = permutest(emp_clus,perm_clus,true,cluster_p,cluster_nperm,false,cluster_n);

% Transpose depending on output orientation
 if size(c_TGM{1},2) > size(c_TGM{1},1)
    p_TGM = p_TGM';
    c_TGM = c_TGM';
    for i = 1:numel(c_TGM)
        c_TGM{i} = c_TGM{i}';
    end
 end

% Find significant clusters
sig_clus = find(p_TGM<=cluster_p);
nsig_clus = find(p_TGM>cluster_p);

sig_TGM = c_TGM(sig_clus); % Filter significant clusters
nsig_TGM = c_TGM(nsig_clus); % Filter significant clusters

% Check whether there are any clusters, and if those clusters are
% significant.
if isempty(sig_clus) && isempty(nsig_clus)~=1
    disp('No significant clusters detected in the TGM');
    clust_sig_ind = [];
elseif isempty(sig_clus) && isempty(nsig_clus)
    disp('No clusters detected in the TGM');
    clust_sig_ind = [];
    clust_nsig_ind = [];
elseif isempty(sig_clus)~=1 && isempty(nsig_clus)==1
    clust_nsig_ind = [];
end

% Vectorize significant clusters
for i = 1:numel(sig_TGM)
    if i == 1
        clust_sig_ind = sig_TGM{i};
    else
        clust_sig_ind = [clust_sig_ind;sig_TGM{i}];
    end
end

% Vectorize non-significant clusters
for i = 1:numel(nsig_TGM)
    if i == 1
        clust_nsig_ind = nsig_TGM{i};
    else
        clust_nsig_ind = [clust_nsig_ind;nsig_TGM{i}];
    end
end

%% PLOTTING
% Plot Diagonal option
if strcmp(mapmethod,'diag')
    % Calculate average diagonal
    diagavg = diag(TGMavg);
    diagavg_Z = squeeze(mean(emp_clus,2));
    
    figure;subplot(5,5,1:10);hold on;
    bt_figure('halflong');
    plot(diagavg_Z,'LineWidth',3,'Color',bt_colorscheme('diagonal'));
    
    try % For old Matlab versions
        yline(0,'LineWidth',1.5,'Color',[0.6 0.6 0.6]);
    catch
        vline(0);
    end
    
    xlim([0,numel(diagavg_Z)]);
    xlabel('Test data (bin)')
    ylabel('Z-value');
    
    if cluster_smooth == 0
        title('Z-scored empirical TGM Diagonal');
    else
        title('Smoothed Z-scored empirical TGM Diagonal');
    end
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    
    % Plot points a little below lowest Z
    lim = min(diagavg_Z(:))*1.05;
    
    % Plot significant and non-significant datapoints
    if isempty(clust_sig_ind) && isempty(clust_nsig_ind)~=1
        c1 = plot(clust_nsig_ind,lim,'o','linewidth', 0.1, 'color',bt_colorscheme('nsigcluster'),'MarkerFaceColor',bt_colorscheme('nsigcluster'));hold on;
        l1 = legend(c1(1),'Non-significant clusters');
    elseif isempty(clust_sig_ind)~=1 && isempty(clust_nsig_ind)
        c1 = plot(clust_sig_ind,lim,'o','linewidth', 0.1, 'color',bt_colorscheme('sigcluster'),'MarkerFaceColor',bt_colorscheme('sigcluster'));hold on;
        l1 = legend(c1(1),'Significant clusters');
    elseif isempty(clust_sig_ind)~=1 && isempty(clust_nsig_ind)~=1
        c1 = plot(clust_sig_ind,lim,'o','linewidth', 0.1, 'color', bt_colorscheme('sigcluster'),'MarkerFaceColor',bt_colorscheme('sigcluster'));hold on;
        c2 = plot(clust_nsig_ind,lim,'o','linewidth', 0.1, 'color', bt_colorscheme('nsigcluster'),'MarkerFaceColor',bt_colorscheme('nsigcluster'));hold on;
        l1 = legend([c1(1),c2(1)],'Significant clusters', 'Non-significant clusters');
    end
    
    % Plot regular average diagonal;
    subplot(5,5,20:25);hold on;
    plot(diagavg,'LineWidth',3,'Color',bt_colorscheme('diagonal'));
    
    try % For old Matlab versions
        yline(0.5,'LineWidth',1.5,'Color',[0.6 0.6 0.6]);
    catch
        vline(0.5);
    end
    
    xlim([0,numel(diagavg_Z)]);
    xlabel('Test data (bin)')
    ylabel(MVPAcfg.metric);
    title('average empirical TGM Diagonal (classifier time course)');
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    l1.FontSize = bt_plotparams('FontSizeLegend');
    
    % Plot TGM option
else
    % Get average z-scored TGM
    TGMavg_Z = squeeze(mean(emp_clus,3));
    
    % Pre-allocate TGM mask and p-value dist
    TGM_sig_mask = zeros(size(TGMavg,1),size(TGMavg,2));
    TGM_nsig_mask = zeros(size(TGMavg,1),size(TGMavg,2));
    
    % Enter significant clusters into mask
    TGM_sig_mask(clust_sig_ind) = true;
    TGM_nsig_mask(clust_nsig_ind) = true;
    
    % Plot clusters in mean empirical TGM
    figure;hold on;
    bt_figure('cluster');
    subplot(10,2,[3 5 7 9 11 13 15]);
    cfg_plot = [];
    cfg_plot.colorbar = 1;
    cfg_plot.colorbar_location = 'EastOutside';
    cfg_plot.colorbar_title = 'Z-value';
    cfg_plot.climzero = 0;
    mv_plot_2D(cfg_plot,TGMavg_Z);hold on;
    colormap(bt_colorscheme('TGM'));freezeColors;
    if cluster_smooth == 0
        title('Z-scored average TGM (clusters)')
    else
        title('Smoothed Z-scored average TGM (clusters)')
    end
    xlabel('Test data (bin)')
    ylabel('Training data (bin)')
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    
    contour(1:size(TGMavg,1),1:size(TGMavg,2),TGM_sig_mask,1,'linewidth', 2, 'color',bt_colorscheme('sigcluster'));hold on;
    contour(1:size(TGMavg,1),1:size(TGMavg,2),TGM_nsig_mask,1,'linewidth', 2, 'color',bt_colorscheme('nsigcluster'));
    l1 = legend('Significant clusters', 'Non-significant clusters');
    
    subplot(10,2,[4 6 8 10 12 14 16]);
    cfg_plot.colorbar_title = MVPAcfg.metric;
    cfg_plot.climzero = 0.5;
    mv_plot_2D(cfg_plot,TGMavg);hold on;
    colormap(bt_colorscheme('TGM'));freezeColors;
    title('Average empirical TGM')
    xlabel('Test data (bin)')
    ylabel('Training data (bin)')
    % Adapt font
    set(gca,'FontName',bt_plotparams('FontName'));
    set(gca,'FontSize',bt_plotparams('FontSize'));
    l1.FontSize = bt_plotparams('FontSizeLegend');
end

%% Adapt output variable
% Add cluster correction information to output
stats2.clusters.clusters_inds = c_TGM;
stats2.clusters.clusters_p = p_TGM;
end