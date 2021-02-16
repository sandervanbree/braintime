function [stats2] = bt_TGMstatslevel2(config, stats1)
% Performs second level statistics of recurrence power spectra, and tests
% which TGM datapoints are significantly higher than the distribution of
% permuted TGMs.
%
% - Second level statistics uses the first level permutated data. Numperm2
% times, grab a random power spectrum of every participants' numperm1 pool
% of permutated data and average. Then, for every frequency, the pvalue is
% the percentile of the average empirical data's recurrence power spectra
% in the null distribution of the permuted data's recurrence power spectra.
%
% - TGM significance testing compares each empirical TGM datapoint to
% the same datapoint in the permuted distribution and tests whether
% classifier performance exceeds 5%. Empirical and permuted  TGMs are
% normalized using the mean and standard deviation of the permuted TGMs.
%
% Use:
% [pval] = bt_TGMstatslevel2(config, stats1), where stats1 has one cell
% for each participant.
%
% Input Arguments:
% config
%   - numperms2      % Number of permutations on the second statistical
%                    % level. Default: 1000000.
%                    %
%   - multiplecorr   % Correction method for multiple comparisons.
%                    % 'none' (default), 'fdr' for adjusted p-values based
%                    % on the False Discovery Rate, and 'bonferroni' for
%                    % Bonferroni correction.
%                    %
% stats1             % Output structure obtained using bt_TGMstatslevel1.
%                    % Each cell contains one participants' empirical TGM,
%                    % permuted TGMs, recurrence power spectrum of the
%                    % empirical and permuted TGM, and the associated
%                    % frequency vector.
%                    %
% Output:            %
% pval               % Three row vectors. The first row vector contains the
%                    % pvalues of the statistical comparison between the
%                    % empirical power spectrum (averaged across
%                    % participants) and the numperm2 permuted datasets.
%                    % The second row vector contains the tested
%                    % frequencies.

%% Get information
numsubj = numel(stats1);                             % Number of participants
numperms1 = size(stats1{1}.permspec,1);             % Number of first level permutations
if isfield(config,'numperms2')                       % Number of second level permutations
    numperms2 = config.numperms2;
else
    numperms2 = 100000;                              % Default to 100000
end
refdimension = stats1{1}.refdimension;               % Find out the reference dimension (clock or brain time)

% Sanity check
if numperms2 <100000
    numperms2 = 100000;
    warning(['Due to high variance in this statistical procedure, the number of second level permutations has been increased from ', num2str(numperms2),' to 100000.']);
end

%% STEP 1: SECOND LEVEL STATISTICS OF RECURRENCE POWER SPECTRA
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

%% Get average recurrence power spectrum across participants
% PS_emp = zeros(numsubj,numel(f));
for subj = 1:numsubj
    PS_emp(subj,:) = stats1{subj}.empspec;
end
PS_emp_avg = mean(PS_emp,1);

%% Second level statistics
% Pre-allocate
perm1PS = zeros(numsubj,mode(nfbins));
perm2PS = zeros(numperms2,mode(nfbins));
progbar = 0:1000:numperms2;

% Loop through numperms2
for perm2=1:numperms2
    if ismember(perm2,progbar) % Print progress
        disp(strcat((num2str(round((perm2/numperms2)*100,2))),'% of second level permutations completed (for second-level stats)'));
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
figure;hold on;set(gcf, 'WindowState', 'maximized'); % create full screen figure

if strcmp(refdimension.dim,'braintime') % Only for brain time, separate warped frequency
    subplot(1,10,1:3)
    
    % Plot the warped freq
    wfreq_emp = PS_emp_avg(:,wfreq_i); % Empirical data at warped freq
    wfreq_perm = perm2PS(:,wfreq_i); % Permuted data at warped freq
    plot(2,wfreq_emp,'o','MarkerSize',6,'MarkerEdgeColor','blue','MarkerFaceColor','blue');hold on; %Plot marker of empirical power
    
    if exist('wfreq_half_i','var') == 1 %Plot half the warped freq?
        half_wfreq_emp = PS_emp_avg(:,wfreq_half_i);
        half_wfreq_perm = perm2PS(:,wfreq_half_i);
        plot(1,half_wfreq_emp,'o','MarkerSize',6,'MarkerEdgeColor','blue','MarkerFaceColor','blue');hold on; %Plot marker of empirical power
    else
        half_wfreq_emp = zeros(size(wfreq_emp));
        half_wfreq_perm = zeros(size(wfreq_perm));
    end
    
    if exist('wfreq_double_i','var') == 1 %Plot double the warped freq?
        double_wfreq_emp = PS_emp_avg(:,wfreq_double_i);
        double_wfreq_perm = perm2PS(:,wfreq_double_i);
        plot(3,double_wfreq_emp,'o','MarkerSize',6,'MarkerEdgeColor','blue','MarkerFaceColor','blue');hold on; %Plot marker of empirical power
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
        l2 = legend(h([numel(h) 3]),'Empirical (emp) recurrence power','Permuted (perm) recurrence power');
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
        l2 = legend([h(2),toMark(1)],'Empirical (emp) recurrence power','Permuted (perm) recurrence power');
        set(l2,'Location','best');
    end
    
    % Set up axes
    ylabel('Recurrence power');
    title('2nd level recurrence')
    
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
    set(gca,'FontName','Arial')
    set(gca,'FontSize',15)

    title('2nd level recurrence')
    
    % Adapt font
    set(gca,'FontName','Arial')
    set(gca,'FontSize',16)
    
    % Prepare new subplot
    subplot(1,10,5:10)
    hold on
end

yyaxis left
p1 = plot(f,PS_emp_avg,'LineStyle','-','LineWidth',3,'Color','b','Marker','none'); %Mean across 2nd level perms
p2 = plot(f,perms2PS_avg,'LineStyle','-','LineWidth',2,'Color',[0.3 0.3 0.3],'Marker','none'); %Mean across 2nd level perms
xlabel('Recurrence frequency')
ylabel('Mean power across participants')

% Plot confidence interval
c2 = plot(f,low_CI,'LineStyle','-','LineWidth',0.5,'Color','k','Marker','none');
c3 = plot(f,hi_CI,'LineStyle','-','LineWidth',0.5,'Color','k','Marker','none');
p3 = patch([f fliplr(f)],[low_CI' fliplr(hi_CI')], 1,'FaceColor', 'black', 'EdgeColor', 'none', 'FaceAlpha', 0.15);

% p-value axis
yyaxis right
p4 = plot(f,logpval,'LineStyle','-','LineWidth',2,'Color',[0.7 0.2 0.2]);
p4.Color(4) = 0.3;
ylim([-0.05,4])
ylabel('-log10 p-value')

% plot star at every significant frequency
yyaxis left
sigind = find(pval_corr<=0.05);
p5 = plot(f(sigind),0,'r*','MarkerSize',10,'LineWidth',1.5);

% legend
l2 = legend([p1 p2 p3 p4 p5(1)],{'Average emp spectrum','Average perm spectrum', 'Conf. interv. perm spectrum', '-log10 p-value', 'p<=0.05'});
set(l2,'Location','best')

% add title
title('Recurrence power spectra (2nd level stats)');

% Adapt font
set(gca,'FontName','Arial')
set(gca,'FontSize',16)

%% Adapt output variable
% Add freq and pval information to output
stats2.pval_corr = pval_corr;
stats2.frequency = f;

% Print results for brain time
disp('See the output structure for:')
disp('(1) the (corrected) p-values per frequency');
disp('(2) the range of tested frequencies');
disp('(3) for brain time, the (corrected) p-value of 0.5x, 1x, and 2x the warped frequency');

% %% STEP 2: TEST SIGNIFICANCE OF TGM
% % Pre-allocate
% empset_TGM = zeros(numel(stats1),size(stats1{1}.empTGM,1),size(stats1{1}.empTGM,2));
% permset_TGM = zeros(numel(stats1),size(stats1{1}.permTGM,1),size(stats1{1}.permTGM,2),size(stats1{1}.permTGM,3));
% 
% % Put each participant's empirical and permuted data into one matrix, whilst z-scoring based on permuted data
% disp('To test for significant points in the TGM, the toolbox normalizes to');
% disp('the mean and standard deviation of permuted TGMs');
% for i = 1:numel(stats1)
%     % Calculate mean and std of the permuted TGMs
%     mn = mean(stats1{i}.permTGM(:));
%     sd = std(stats1{i}.permTGM(:));
%     
%     % Subtract the mean and divide by std to normalize across participants
%     empset_TGM(i,:,:) = (stats1{i}.empTGM-mn)./sd;
%     permset_TGM(i,:,:,:) = (stats1{i}.permTGM-mn)./sd;
% end
% 
% % Create average empirical TGM
% avgemp = squeeze(mean(empset_TGM,1));
% 
% % Pre-allocate
% tempperm = zeros(numsubj,size(permset_TGM,3),size(permset_TGM,4));
% distX = zeros(numperms2,size(permset_TGM,3),size(permset_TGM,4));
% 
% % Loop through numperms2
% for perm2 = 1:numperms2
%     if ismember(perm2,progbar) % Print progress
%         disp(strcat((num2str(round((perm2/numperms2)*100,2))),'% of second level permutations completed (for TGM significance testing)'));
%     end
%     for s = 1:numsubj      % For each participant
%         idx = randperm(numperms1,1); % Grab a random first permutation
%         tempperm(s,:,:) = squeeze(permset_TGM(s,idx,:,:));
%     end
%     distX(perm2,:,:) = squeeze(mean(tempperm,1));
% end
% 
% nrow = size(distX(1,:,:),2);
% ncol = size(distX(1,:,:),3);
% 
% % Get p-values of each TGM datapoint
% pval=zeros(nrow,ncol);
% for r=1:nrow
%     for c=1:ncol
%         pval(r,c)=numel(find(distX(:,r,c)>=avgemp(r,c)))/numperms2;
%     end
% end
% 
% % Multiple comparison correction
% pval_vec=pval(:);
% [~,~,~,pval_corr] = fdr_bh(pval_vec,0.05,'dep');
% pval_corr = reshape(pval_corr,nrow,ncol);
% 
% % Logical mask based on significant points
% sigmask=zeros(nrow,ncol);
% for r=1:nrow
%     for c=1:ncol
%         if pval_corr(r,c)<=0.05 % Only let through significant points
%             sigmask(r,c)=1;
%         end
%     end
% end
% 
% sig_i = find(sigmask==1); % Find only significant points
% maskedTGM = nan(nrow,ncol); % Prepare NaN TGM
% maskedTGM(sig_i) = avgemp(sig_i); % Replace by values of average emp TGM
% 
% figure;
% pcolor(1:nrow,1:ncol,maskedTGM);
% shading interp;
% xlabel('Data points (train)');
% ylabel('Data points (test)');
% title('TGM data points significantly higher than permuted (FDR-corrected)');
