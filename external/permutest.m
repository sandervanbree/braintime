function [clusters, p_values, t_sums, permutation_distribution ] = permutest( trial_group_1, trial_group_2, dependent_samples, ...
    p_threshold, num_permutations, two_sided, num_clusters )
% Permutation test for dependent or independent measures of 1-D or 2-D data. 
% Based on Maris & Oostenveld 2007 for 1-D and 2-D vectors. The test 
% statistic is T-Sum - the total of t-values within a cluster of contingent
% above-threshold data points. 
% See: Maris, E., & Oostenveld, R. (2007). Nonparametric statistical 
% testing of EEG-and MEG-data. Journal of Neuroscience Methods, 164(1), 
% 177?190. https://doi.org/10.1016/j.jneumeth.2007.03.024
% 
% Important notes: 
% * Make sure you understand whether you should be using a test of
% dependent or independent samples (a within or between subjects test,
% respectively). Also make sure you know if you want a one-sided or
% two-sided test. 
% * The number of permutation sets a minimal boundary for the resulting
% p-value. This boundary is 1/(nP+1). If a very low p-value is needed, e.g.
% to survive a multiple-comparisons correction, use a sufficiently high
% number of permutations (this is also necessary in order to estimate the
% p-value with sufficient accuracy unless there is a very small number of
% trials). 
% * When runnning a two-sided test, there is no need to retroactively
% correct the p-value, as this more lenient assumption is already reflected 
% in the way the null-hypothesis distribution is constructed. 
% 
% Syntax: 
% [clusters, p_values, t_sums, permutation_distribution ] = PERMUTEST(x,y) 
% runs a cluster-based permutation based for dependent samples, testing for 
% differences between x and y. x and y should be 2D or 3D matrices, where 
% the last dimension is the dimension across which the test variance is 
% defined (that is, trials or subjects). clusters is a cell array with each 
% cell holding the vector/matrix indexes corresponding to a specific 
% cluster (sorted by magnitude). p_values is the permutation p-value 
% corresponding to each cluster. t_sums is the sum of t-values of all data 
% points comprising each cluster. distribution is the T-Sum permutation 
% distribution (note that it will include many zero values; these 
% correspond to all the permutations where no above-threshold clusters were 
% found). 
% PERMUTEST(x,y,d), where d is set to true or false, determines whether the
% test is for dependent samples. If set to false, it is assumed that x and
% y are independent (non-paired) data sets. In this case, x and y can have
% a different number of trials (though all other dimensions should be equal
% in size). Default value is true. 
% PERMUTEST(x,y,d,p) sets p as the p-value threshold
% below which a data point is part of a cluster. Lower values mean that the 
% test is sensitive to narrow, stronger effects; higher values make the 
% test sensitive to broad, weaker effects.  The p-value is translated to a 
% t-value for the purpose of this test. 
% PERMUTEST(x,y,d,p,nP) sets nP as the number of permutations. The function
% will check whether this number can be supported for the given number of
% trials and test type. 
% PERMUTEST(x,y,d,p,nP,t), where t is set to true or false, determines
% whther the test is a two-sided test. If set to true, negative differences
% will be detected as well. Default value is false (one-sided test). 
% PERMUTEST(x,y,d,p,nP,t,nC) sets nC as the maximal number of significant
% clusters to be detected. Default value is inf (all existing clusters will
% be tested against the H0 distribution). 
%
% Written by Edden M. Gerber, lab of Leon Y. Deouell, 2014
% Send bug reports and requests to edden.gerber@gmail.com
% INITIALIZE
% Set optional arguments:
if nargin < 7 || isempty(num_clusters)
    num_clusters = inf;
end
if nargin < 6 || isempty(two_sided)
    two_sided = false;
end
if nargin < 5 || isempty(num_permutations)
    num_permutations = 10^4;
end
if nargin < 4 || isempty(p_threshold)
    p_threshold = 0.05;
end
if nargin < 3 || isempty(dependent_samples)
    dependent_samples = true;
end
if nargin < 2
    error('Not enough input arguments');
end
% Check input dimensions:
if ismatrix(trial_group_1)
    [num_data_points_1_x, num_trials_1] = size(trial_group_1);
    num_data_points_1_y = 1;
elseif ndims(trial_group_1) == 3
    [num_data_points_1_x, num_data_points_1_y, num_trials_1] = size(trial_group_1);
else
    error('Input variable "trial_group_1" needs to be 2D or 3D');
end
if ismatrix(trial_group_2)
    [num_data_points_2_x, num_trials_2] = size(trial_group_2);
    num_data_points_2_y = 1;
elseif ndims(trial_group_2) == 3
    [num_data_points_2_x, num_data_points_2_y, num_trials_2] = size(trial_group_2);
else
    error('Input variable "trial_group_2" needs to be 2D or 3D');
end
if dependent_samples
    num_trials = num_trials_1;
    if ~isequal(num_data_points_1_x,num_data_points_2_x) || ~isequal(num_data_points_1_y,num_data_points_2_y) ...
            || ~isequal(num_trials_1, num_trials_2)
        error('Size of all dimensions should be identical for two dependent samples');
    end
else
    if ~isequal(num_data_points_1_x,num_data_points_2_x) || ~isequal(num_data_points_1_y,num_data_points_2_y)
        error('Size of all dimensions but the last one (corresponding to the number of trials) should be identical for two independent samples');
    end
end
num_data_points_x = num_data_points_1_x;
num_data_points_y = num_data_points_1_y;
% Check that number of requested permutations is possible:
if dependent_samples
    % In each permutation, each paired sample could be flipped, meaning
    % that there are 2^num_trials possible permutation
    max_num_permutations = 2^num_trials;
    if num_permutations > max_num_permutations
        warning('With %d paired trials, only %d permutations are possible. Using this value instead of %d.',...
            num_trials, max_num_permutations, num_permutations);
        num_permutations = max_num_permutations;
    end
else
    % For independent samples, the number of possible permutations is
    % nchoosek(num_trials_1+num_trials_2, num_trials_1) - since it is
    % assumed that the same trial group sizes are maintained across all
    % permutations. 
    warning('off','MATLAB:nchoosek:LargeCoefficient'); % no need to alarm anyone with this warning
    max_num_permutations = nchoosek(num_trials_1+num_trials_2, num_trials_1);
    warning('on','MATLAB:nchoosek:LargeCoefficient')
    if num_permutations > max_num_permutations
        warning('With %d trials in group1 and %d trials in group2, only %d permutations are possible. Using this value instead of %d',...
            num_trials_1, num_trials_2, max_num_permutations, num_permutations);
        num_permutations = max_num_permutations;
    end
end
% Initialize output variables
clusters = cell(1);
p_values = ones(1);
% Compute t-value threshold from p-value threshold
if dependent_samples
    tThreshold = abs(tinv(p_threshold, num_trials-1));
else
    tThreshold = abs(tinv(p_threshold, num_trials_1+num_trials_2-1));
end
% PRODUCE PERMUTATION VECTORS
if num_permutations < max_num_permutations / 1000
    % there are at least 1000 times more possible permutations than the 
    % requested number, so draw permutation patterns randomly.                             
    if dependent_samples
        % Dependent samples: randomly generate vectors of 1's and -1's. In
        % each permutation, each pair's difference will be multipled by
        % either 1 or -1. 
        permutation_vectors = round(rand(num_trials,num_permutations)) * 2 - 1; 
    else
        % Independent samples: randomly generate vectors of 1's and 2's. In
        % each pemutation, this vector will dictate to which group each
        % trial belongs (the number of trials from each group is identical
        % to that of the original samples). 
        permutation_vectors = ones(num_trials_1+num_trials_2, num_permutations);
        for p = 1:num_permutations
            idx = randperm(num_trials_1+num_trials_2, num_trials_2);
            permutation_vectors(idx, p) = 2;
        end
    end
else
    % The number of requested permutations is close to the maximum which 
    % can be produced by the number of trials, so permutation sequences are 
    % not randomly drawn but generated explicitly without repetition. 
    if dependent_samples
        % Dependent samples: randomly generate vectors of 1's and -1's. In
        % each permutation, each pair's difference will be multipled by
        % either 1 or -1. 
        
        % Initialize variable:
        permutation_vectors = NaN(num_trials,num_permutations);
        % generate a non-repeating list of permutations by taking a list of
        % non-repating numbers (up to nPermutations), and translating them 
        % into binary (0's and 1's): 
        rndB = dec2bin(randperm(2^num_trials,num_permutations)-1);
        % if the leading bit of all the numbers is 0, it will be truncated, so
        % we need to fill it in:
        nBits = size(rndB,2);
        if nBits < num_trials
            rndB(:,(num_trials-nBits+1):num_trials) = rndB;
            rndB(:,1:num_trials-nBits) = '0';
        end
        % translate the bits into -1 and 1:
        for ii=1:numel(permutation_vectors)
            permutation_vectors(ii) = str2double(rndB(ii)) * 2 - 1;
        end
    else
        % Independent samples: randomly generate vectors of 1's and 2's. In
        % each pemutation, this vector will dictate to which group each
        % trial belongs (the number of trials from each group is identical
        % to that of the original samples). 
        
        % Initialize variable:
        permutation_vectors = ones(num_trials_1+num_trials_2,num_permutations);
        % Generate all possible combinations of num_trials_2 out of
        % (num_trials_1+num_trials_2) entries
        idx_matrix = nchoosek(1:(num_trials_1+num_trials_2),num_trials_2);
        % Randomly draw num_permutations out of them
        idx_matrix = idx_matrix(randperm(size(idx_matrix,1),num_permutations),:);
        % Generate the vectors of 1's and 2's using the above
        for p=1:num_permutations
            permutation_vectors(idx_matrix(p,:),p) = 2;
        end
    end
end
% RUN PRIMARY PERMUTATION
t_value_vector = zeros(num_data_points_x,num_data_points_y);
if dependent_samples
    % Dependent samples: one-sample t-test of the difference between the 
    % two samples (compared to zero)
    if num_data_points_y == 1 
        % one-dimensional trial data
        t_value_vector = simpleTTest(trial_group_1'-trial_group_2',0);
    else
        % two-dimensional trial data
        for ii = 1:num_data_points_y
            t_value_vector(:,ii) = simpleTTest(squeeze(trial_group_1(:,ii,:)-trial_group_2(:,ii,:))',0);
        end
    end
else
    % Independent samples: two-sample t-test of the difference between the
    % 2 sample means
    if num_data_points_y == 1 
        % one-dimensional trial data
        t_value_vector = simpleTTest2(trial_group_1',trial_group_2');
    else
        % two-dimensional trial data
        for ii = 1:num_data_points_y
            t_value_vector(:,ii) = simpleTTest2(squeeze(trial_group_1(:,ii,:))',squeeze(trial_group_2(:,ii,:))');
        end
    end
end
% Find the above-threshold clusters: 
CC = bwconncomp(t_value_vector > tThreshold,4);
cMapPrimary = zeros(size(t_value_vector));
tSumPrimary = zeros(CC.NumObjects,1);
for i=1:CC.NumObjects
    cMapPrimary(CC.PixelIdxList{i}) = i;
    tSumPrimary(i) = sum(t_value_vector(CC.PixelIdxList{i}));
end
if two_sided % Also look for negative clusters
    n = CC.NumObjects;
    CC = bwconncomp(t_value_vector < -tThreshold,4);
    for i=1:CC.NumObjects
        cMapPrimary(CC.PixelIdxList{i}) = n+i;
        tSumPrimary(n+i) = sum(t_value_vector(CC.PixelIdxList{i}));
    end
end
% Sort clusters:
[~,tSumIdx] = sort(abs(tSumPrimary),'descend');
tSumPrimary = tSumPrimary(tSumIdx);
% RUN PERMUTATIONS
% We'll need this for creating shuffled trials for independent samples:
if ~dependent_samples
    all_trials = cat(ndims(trial_group_1), trial_group_1, trial_group_2);
end
permutation_distribution = zeros(num_permutations,1);
for p = 1:num_permutations
    if dependent_samples
        % Dependent samples: use a one-sample t-test of the difference 
        % between the two samples (compared to zero)
        if num_data_points_y == 1 
            % one-dimensional trial data
            D = bsxfun(@times,trial_group_1'-trial_group_2',permutation_vectors(:,p));
            t_value_vector = simpleTTest(D,0);
        else
            % two-dimensional trial data
            for ii = 1:num_data_points_y
                D = squeeze(trial_group_1(:,ii,:)-trial_group_2(:,ii,:))';
                D = bsxfun(@times,D,permutation_vectors(:,p));
                t_value_vector(:,ii) = simpleTTest(D,0);
            end
        end
    else
        % Independent samples: use a two-sample t-test of the difference 
        % between the 2 sample means
        if num_data_points_y == 1 
            % one-dimensional trial data    
            p_trial_group_1 = all_trials(:,permutation_vectors(:,p)==1);
            p_trial_group_2 = all_trials(:,permutation_vectors(:,p)==2);
            t_value_vector = simpleTTest2(p_trial_group_1',p_trial_group_2');
        else
            % two-dimensional trial data
            p_trial_group_1 = all_trials(:,:,permutation_vectors(:,p)==1);
            p_trial_group_2 = all_trials(:,:,permutation_vectors(:,p)==2);
            for ii = 1:num_data_points_y
                t_value_vector(:,ii) = simpleTTest2(squeeze(p_trial_group_1(:,ii,:))',squeeze(p_trial_group_2(:,ii,:))');
            end
        end
    end
    
    % Find clusters:
    CC = bwconncomp(t_value_vector > tThreshold,4);
    tSum = zeros(CC.NumObjects,1);
    for i=1:CC.NumObjects
        tSum(i) = sum(t_value_vector(CC.PixelIdxList{i}));
    end
    if two_sided % Also look for negative clusters
        n = CC.NumObjects;
        CC = bwconncomp(t_value_vector < -tThreshold,4);
        for i=1:CC.NumObjects
            tSum(n+i) = sum(t_value_vector(CC.PixelIdxList{i}));
        end
    end
    if isempty(tSum)
        permutation_distribution(p) = 0;
    else
        [~,idx] = max(abs(tSum));
        permutation_distribution(p) = tSum(idx);
    end
end
%% DETERIMNE SIGNIFICANCE
for clustIdx = 1:min(num_clusters,length(tSumPrimary))
    if two_sided
        ii = sum(abs(permutation_distribution) >= abs(tSumPrimary(clustIdx)));
    else
        ii = sum(permutation_distribution >= tSumPrimary(clustIdx));
    end
    
    clusters{clustIdx} = find(cMapPrimary == tSumIdx(clustIdx));
    p_values(clustIdx) = (ii+1) / (num_permutations+1);
end
% return regular arrays if only one cluster is requested
t_sums = tSumPrimary;
if num_clusters == 1
    clusters = clusters{1};
    permutation_distribution = permutation_distribution{1};
end
end
function [t, df] = simpleTTest(x,m)
    %TTEST  Hypothesis test: Compares the sample average to a constant.
    %   [STATS] = TTEST(X,M) performs a T-test to determine
    %   if a sample from a normal distribution (in X) could have mean M.
    %  Modified from ttest function in statistical toolbox of Matlab
    % The  modification is that it returns only t value and df. 
    %  The reason is that calculating the critical value that
    % passes the threshold via the tinv function takes ages in the original
    % function and therefore it slows down functions with many
    % iterations.
    % Written by Leon Deouell. 
    % 
    % Modified by Edden Gerber 19.11.2013: Added support for x being a matrix, where columns are
    % observations and rows are variables (output is a vector). 
    if nargin < 1, 
        error('Requires at least one input argument.'); 
    end
    if nargin < 2
        m = 0;
    end
    samplesize  = size(x,1);
    xmean = sum(x)/samplesize; % works faster then mean
    % compute std  (based on the std function, but without unnecessary stages
    % which make that function general, but slow (especially using repmat)
    xc = bsxfun(@minus,x,xmean);  % Remove mean
    xstd = sqrt(sum(conj(xc).*xc,1)/(samplesize-1));
    ser = xstd ./ sqrt(samplesize);
    tval = (xmean - m) ./ ser;
    % stats = struct('tstat', tval, 'df', samplesize-1);
    t = tval;
    df = samplesize-1;
end
function [t, df] = simpleTTest2(x1,x2)
% A 2-sample t-test which computes only the t-value (and degrees of
% freedom), skipping the very time-consuming p-value computation. This
% function is good for permutation tests which need to compute t-values a
% large number of times as fast as possible. 
% If using input matrices, different variables should be along the rows
% (each column is a single variable).
%
% Written by Edden Gerber, March 2016. 
if nargin < 2, 
    error('Requires at least two input argument.'); 
end
n1  = size(x1,1);
n2  = size(x2,1);
xmean1 = sum(x1)/n1; % works faster then mean
xmean2 = sum(x2)/n2; % works faster then mean
% compute std (based on the std function, but without unnecessary stages
% which make that function general, but slow (especially using repmat)
xc = bsxfun(@minus,x1,xmean1);  % Remove mean
xstd1 = sqrt(sum(conj(xc).*xc,1)/(n1-1));
xc = bsxfun(@minus,x2,xmean2);  % Remove mean
xstd2 = sqrt(sum(conj(xc).*xc,1)/(n2-1));
sx1x2 = sqrt(xstd1.^2/n1 + xstd2.^2/n2);
t = (xmean1 - xmean2) ./ sx1x2;
df = n1+n2-2;
end