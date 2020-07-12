function [pval] = bt_TGMstatslevel2(config, stats1)
% Acquire second level statistics using the first level permutated
% data. Numperm2 times, grab a random power spectrum of every participants'
% numperm1 pool of permutated data and average. Then, for every frequency,
% the pvalue is the percentile of the average empirical data's power in the
% null distribution of the permuted data's power.
%
% Use:
% [pval] = bt_TGMstatslevel2(config, stats1), where stats1 has one cell
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
% stats1             % Output structure obtained using bt_TGMstatslevel1.
%                    % Each cell contains one participants' power spectra
%                    % the empirical data, permutation data, the associated
%                    % frequency vector, and the power at the mode freq.
%                    % 
% Output:            % 
% pval               % Two row vectors. The first row vector contains the
%                    % pvalues of the statistical comparison between the
%                    % empirical power spectrum (averaged across
%                    % participants) and the numperm2 shuffled datasets.
%                    % The second row vector contains the tested
%                    % frequencies.

%% Get information
numsubj = numel(stats1); %number of participants
numperms1 = size(stats1{1}.shuffspec,1); %number of first level permutations
numperms2 = config.numperms2; %number of second level permutations
f = stats1{1}.f; %frequency vector
% fullspec_emp = stats1.empspec; %full power spectrum of empirical data
% fullspec_shuff = stats1.shuffspec; %full power spectrum of shuffled data
% modepow_emp = stats1.empmode; %mode power of empirical data
% modepow_shuff = stats1.shuffmode; %mode power of shuffled data

%% Get average recurrence power spectrum across participants
PS_emp = zeros(numsubj,numel(f));
for subj = 1:numsubj
    PS_emp(subj,:) = stats1{subj}.empspec;
end
PS_emp = mean(PS_emp,1);

%% Second level statistics
for perm2=1:numperms2
    fprintf('Second level permutation number %i\n', perm2);
    for subj = 1:numsubj
        idx=randperm(numperms1,1); %randomly grab with replacement
        perm1PS(subj,:) = stats1{subj}.shuffspec(idx,:);        
    end
    perm2PS(perm2,:) = mean(perm1PS,1);
end

% Mean across all 2nd level permutations (for plotting)
perms2PS_avg = mean(perm2PS,1);

% Maximum difference between emp and shuff
[~, max_f] = max(PS_emp-perms2PS_avg);

%% Calculate p-values
for f_ind = 1:numel(f)
pval(f_ind) = numel(find(perm2PS(:,f_ind)>=PS_emp(f_ind)))/numperms2;
end

%% Multiple comparisons correction
if isfield(config,'multiplecorr')
    if strcmp(config.multiplecorr,'fdr')
        [~,~,~,pval] = fdr_bh(pval);   
    elseif strcmp(config.multiplecorr,'bonferroni')
        pval = pval*numel(pval);
    end
end

logpval = -log10(pval);
logpval(isinf(logpval)) = 4; %cap logpval on 4.

%% Plot results
% relationship empirical and shuffled amp at all freq
figure; hold on
yyaxis left
p1 = plot(f,PS_emp,'LineStyle','-','LineWidth',3,'Color','b'); %Mean across 2nd level perms
p2 = plot(f,perms2PS_avg,'LineStyle','-','LineWidth',2,'Color',[0.3 0.3 0.3]); %Mean across 2nd level perms
p3 = line([f(max_f) f(max_f)], [0 max(PS_emp)],'color',[1 0 1],'LineWidth',3); %Line at warped freq
p3.Color(4) = 0.45;
xlabel('Recurrence frequency')
ylabel('Mean power across participants')

% p-value axis
yyaxis right
p4 = plot(f,logpval,'LineStyle','-','LineWidth',2,'Color',[0.7 0.2 0.2]);
p4.Color(4) = 0.25;
ylabel('-log10 p-value')

% plot star at every significant frequency
yyaxis left
sigind = find(pval<=0.05);
p5 = plot(f(sigind),PS_emp(sigind),'r*','MarkerSize',10,'LineWidth',1.5);

% legend
legend([p1 p2 p3 p4 p5],{'Average emp spectrum','Average perm spectrum','Largest diff emp and perm spectrum','-log10 pvalue', 'p<=0.05'});

%% Adapt output variable
% Add freq information to pval vector
pval(2,:) = f;
