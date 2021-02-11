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
PS_emp = mean(PS_emp,1);

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
    pval(f_ind) = numel(find(perm2PS(:,f_ind)>=PS_emp(f_ind)))/numperms2;
    
    % Confidence interval
    low_CI(f_ind,:) = prctile(perm2PS(:,f_ind),2.5);
    hi_CI(f_ind,:) = prctile(perm2PS(:,f_ind),97.5);
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

% For brain time referenced data, exempt the warped frequency from multiple testing correction
if strcmp(refdimension.dim,'braintime')
wfreq = nearest(f,1); %Find the warped frequency (1 Hz)
pval_corr(wfreq) = pval(wfreq);
fprintf(1, '\n');
disp('##############################################################################%########')
disp('The p-value at the warped frequency (1 Hz) is exempted from multiple testing correction');
disp('#######################################################################################')
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
plot(PS_emp(wfreq),'o','MarkerSize',6,'MarkerEdgeColor','blue','MarkerFaceColor','blue'); %Plot marker of empirical power
violinplot(perm2PS(:,wfreq),'test','ShowData',false,'ViolinColor',[0.8 0.8 0.8],'MedianColor',[0 0 0],'BoxColor',[0.5 0.5 0.5],'EdgeColor',[0 0 0],'ViolinAlpha',0.8);

% Set legend
h = get(gca,'Children');
l2 = legend(h([9 3]),'Empirical (emp) recurrence power','Permuted (perm) recurrence power');
set(l2,'Location','best');

% Set up axes
ylabel('Recurrence power');
xticklabels(' ');
xlabel('Warped frequency (1 Hz)');
title(['P-value at warped frequency (1Hz): p = ',num2str(round(pval_corr(wfreq)),6)])

% Adapt font
set(gca,'FontName','Arial')
set(gca,'FontSize',16)

% Prepare new subplot
subplot(1,10,5:10)
hold on
end

yyaxis left
p1 = plot(f,PS_emp,'LineStyle','-','LineWidth',3,'Color','b'); %Mean across 2nd level perms
p2 = plot(f,perms2PS_avg,'LineStyle','-','LineWidth',2,'Color',[0.3 0.3 0.3]); %Mean across 2nd level perms
xlabel('Recurrence frequency')
ylabel('Mean power across participants')

% Plot confidence interval
c2 = plot(f,low_CI,'LineStyle','-','LineWidth',0.5,'Color','k');
c3 = plot(f,hi_CI,'LineStyle','-','LineWidth',0.5,'Color','k');
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
p5 = plot(f(sigind),PS_emp(sigind),'r*','MarkerSize',10,'LineWidth',1.5);

% legend
l2 = legend([p1 p2 p3 p4 p5],{'Average emp spectrum','Average perm spectrum', 'Conf. interv. perm spectrum', '-log10 p-value', 'p<=0.05'});
set(l2,'Location','best')

% add title
title('Recurrence power spectra (2nd level stats)');

% Adapt font
set(gca,'FontName','Arial')
set(gca,'FontSize',16)

%% Adapt output variable
% Add freq information to pval vector
stats2 = [];
stats2.pval = pval;
stats2.pval_corr = pval_corr;
stats2.frequency = f;

%% STEP 2: TEST SIGNIFICANCE OF TGM
% Pre-allocate
empset = zeros(numel(stats1),size(stats1{1}.empTGM,1),size(stats1{1}.empTGM,2));
permset = zeros(numel(stats1),size(stats1{1}.permTGM,1),size(stats1{1}.permTGM,2),size(stats1{1}.permTGM,3));

% Put each participant's empirical and permuted data into one matrix, whilst z-scoring based on permuted data
disp('To test for significant points in the TGM, the toolbox normalizes to');
disp('the mean and standard deviation of permuted TGMs');
for i = 1:numel(stats1)
    % Calculate mean and std of the permuted TGMs
    mn = mean(stats1{i}.permTGM(:));
    sd = std(stats1{i}.permTGM(:));
    
    % Subtract the mean and divide by std to normalize across participants
    empset(i,:,:) = (stats1{i}.empTGM-mn)./sd;
    permset(i,:,:,:) = (stats1{i}.permTGM-mn)./sd;
end    

% Create average empirical TGM
avgemp = squeeze(mean(empset,1));

% Pre-allocate
tempperm = zeros(numsubj,size(permset,3),size(permset,4));
distX = zeros(numperms2,size(permset,3),size(permset,4));

% Loop through numperms2
for perm2 = 1:numperms2
    if ismember(perm2,progbar) % Print progress
        disp(strcat((num2str(round((perm2/numperms2)*100,2))),'% of second level permutations completed (for TGM significance testing)'));
    end
    for s = 1:numsubj      % For each participant
        idx = randperm(numperms1,1); % Grab a random first permutation
        tempperm(s,:,:) = squeeze(permset(s,idx,:,:));
    end
    distX(perm2,:,:) = squeeze(mean(tempperm,1));
end

nrow = size(distX(1,:,:),2); 
ncol = size(distX(1,:,:),3);

% Get p-values of each TGM datapoint
pval=zeros(nrow,ncol);
for r=1:nrow
    for c=1:ncol
        pval(r,c)=numel(find(distX(:,r,c)>=avgemp(r,c)))/numperms2;
    end
end

% Multiple comparison correction
pval_vec=pval(:);
[~,~,~,pval_corr] = fdr_bh(pval_vec,0.05,'dep');
pval_corr = reshape(pval_corr,nrow,ncol);

% Logical mask based on significant points
sigmask=zeros(nrow,ncol);
for r=1:nrow
    for c=1:ncol
        if pval_corr(r,c)<=0.05 % Only let through significant points
            sigmask(r,c)=1;
        end
    end
end

sig_i = find(sigmask==1); % Find only significant points
maskedTGM = nan(nrow,ncol); % Prepare NaN TGM
maskedTGM(sig_i) = avgemp(sig_i); % Replace by values of average emp TGM

figure;
pcolor(1:nrow,1:ncol,maskedTGM);
shading interp;
xlabel('Data points (train)');
ylabel('Data points (test)');
title('TGM data points significantly higher than permuted (FDR-corrected)');
