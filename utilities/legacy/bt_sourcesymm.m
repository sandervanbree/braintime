function [asymmidx,asymmidx_tdist,wavshap] = bt_sourcesymm(wsources,foi,ncycles)
% This is a legacy function that is no longer part of the Brain Time
% Toolbox, but may come in useful.
% Please see tutorial_checksymmetry and bt_checksymmetry for details.
%
% Written by Simon Hanslmayr. Modified by Sander van Bree

ntrls          = size(wsources.trial,2);                  % Fetch number of trials
totsamples     = numel(wsources.time{1,1});               % Total number of samples per trial
nwsources      = size(wsources.trial{1},1);               % Fetch number of warping sources
nfreqs         = numel(foi);                              % Fetch number of frequencies

% Perform bandpass filtering for each frequency
wsources_bp = cell(nfreqs,1);

count = 1;
for centerfreq = foi
    disp(['Filtering warping sources around ',num2str(centerfreq),' Hz']);
    
winwidth = 2;
    
    cfg            = [];
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = [centerfreq-winwidth/2 centerfreq+winwidth/2];
    wsources_bp{count} = ft_preprocessing(cfg,wsources);
    
    count = count+1;
end

%% Prepare basic variables
% Compute basic parameters
deltat   = 1000/wsources.fsample;
tvec     = -ncycles/2:0.01:ncycles/2;                % Set up time vector for average waveshape

% Pre-allocate primary variables
asc_tmp = [];                                        % Ascending flank of cycles (trough to peak)
desc_tmp = [];                                       % Descending flank of cycles (peak to trough)
wavshap_trl = [];                                % The average wave shape for a trial
wavshap = zeros(nfreqs,nwsources,numel(tvec));   % The average wave shape across all trials
asymmidx = zeros(nfreqs,nwsources);              % Asymmetry index from -1 to 1, representing asymmetry direction (ascending or descending flank)
asymmidx_tdist = asymmidx;                       % T-distribution for asymmetry indices, indicating how much they differ from 0.

%% Loop over data
% For every frequency of interest
for f = 1:numel(wsources_bp)
    data_bp_curr = wsources_bp{f};
    
    curr_centerfreq = foi(f);
    % Parameters for this frequency
    win      = round((1000/curr_centerfreq)/deltat)*(ncycles./2); % Window covering 1 cycle
    halfcyc  = ((1000./curr_centerfreq)./2)./deltat;              % Window covering 1/2 cycle
    locwin   = -floor(halfcyc./4):1:ceil(halfcyc./4);             % Window over which local minima and maxima are calculated
    
    % For every warping source
    for w = 1:nwsources
        c = 0;                                       % Counter
        
        disp(['Calculating waveshape and asymmetry metrics for ',num2str(foi(f)),' Hz (warping source ',num2str(w),' out of ',num2str(nwsources),')...']);
        
        % For every trial
        for n=1:ntrls
            trlf   = data_bp_curr.trial{1,n}(w,:);          % Extract current trial (band pass filtered)
            trlr   = wsources.trial{1,n}(w,:);              % Extract current trial (raw)
            
            [~,pk] = findpeaks(trlf(1,:));             % Detect peaks in filtered data
            [~,tr] = findpeaks(trlf(1,:).*-1);         % Detect troughs in filtered data (just flip sign)
            
            % Find local maxima surrounding PEAKS
            for i1 = 1:length(pk)
                if min(pk(i1)+locwin)> 0 && max(pk(i1)+locwin) < totsamples    % As long as we're note before the beginning or after the end
                    locmax     = max(trlr(pk(i1)+locwin));                     % Local maximum value in raw data
                    locmax_i   = find(trlr(pk(i1)+locwin) == locmax);          % Where is the local maximum in the window?
                    pk_ind     = pk(i1)+locwin(locmax_i(1));                   % Correct peak by taking local maximum of raw data
                    L1_wavshap = pk(i1);                                       % Collect peak for average waveshape computation
                    next_tr    = find(tr>pk(i1));                              % What trough succeeds the current peak?
                    
                    % Get succeeding trough from filtered data and extract local minimum
                    if ~isempty(next_tr)                                       % As long as there is another trough
                        next_tr = next_tr(1);
                        if max(tr(next_tr)+locwin) < totsamples
                            locmin    = min(trlr(tr(next_tr)+locwin));
                            locmin_i  = find(trlr(tr(next_tr)+locwin) == locmin);
                            tr_ind    = tr(next_tr)+locwin(locmin_i(1));
                            desc_tmp = [desc_tmp; tr_ind-pk_ind];                 % (t)trough - t(peak) reflects DESCENDING flank
                        end
                    end
                    
                    % Calculate peak-triggered average for wave shape visualization
                    if L1_wavshap-win > 0 && L1_wavshap+win < totsamples
                        c=c+1;
                        wavshap_trl(c,:) = trlr(L1_wavshap-win:L1_wavshap+win); % Extract window from (size winwidth) from raw data
                    end
                end
            end
            
            % Find local minima surrounding TROUGHS - same algorithm, opposite signs
            for i2 = 1:length(tr)
                if min(tr(i2)+locwin)> 0 && max(tr(i2)+locwin) < totsamples    % As long as we're note before the beginning or after the end
                    locmin    = min(trlr(tr(i2)+locwin));                      % Local minimum value in raw data
                    locmin_i  = find(trlr(tr(i2)+locwin)==locmin);             % Where is the local minimum in the window?
                    tr_ind    = tr(i2)+locwin(locmin_i(1));                    % Correct trough by taking local minimum of raw data
                    nextpk    = find(pk>tr(i2));                               % What trough succeeds the current trough?
                    
                    % Get succeeding peak from filtered data and extract local maximum
                    if ~isempty(nextpk)                                        % As long as there is another trough
                        nextpk = nextpk(1);
                        if max(pk(nextpk)+locwin) < totsamples
                            locmax    = max(trlr(pk(nextpk)+locwin));
                            locmax_i  = find(trlr(pk(nextpk)+locwin) == locmax);
                            pk_ind    = pk(nextpk)+locwin(locmax_i(1));
                            asc_tmp   = [asc_tmp; pk_ind-tr_ind];                  % (t)peak - t(trough) reflects ASCENDING flank
                        end
                    end
                end
            end
            % Correct for inconsistencies in number of datapoints for ascending an descending flank
            if numel(asc_tmp)>numel(desc_tmp)
                idx=numel(desc_tmp)+1:numel(asc_tmp);
                asc_tmp(idx)=[];
            elseif numel(desc_tmp)>numel(asc_tmp)
                idx=numel(asc_tmp)+1:numel(desc_tmp);
                desc_tmp(idx)=[];
            end
            
        end
        
        % Average asymmetry index and waveshape across trials warping
        % source and frequency
        currdist = (asc_tmp-desc_tmp)./(halfcyc*2);
        
        % This normalizes the asymmetry index to give a ratio of how much of the
        % cylce is covered by ascending vs descending flank using a value which
        % ranges between -1 and 1.
        asymmidx(f,w) = mean(currdist);
        
        % Get the t-statistic for this frequency and warping source
        [~,~,~,tresult] = ttest(currdist,0,'Tail','both');
        asymmidx_tdist(f,w) = tresult.tstat;
        
        % This calculates the average waveshape and resizes it
        tmp_mn = mean(wavshap_trl,1);
        wavshap(f,w,:) = imresize(tmp_mn,[1 numel(tvec)]);
        
        % Reset variables
        asc_tmp = [];
        desc_tmp = [];
        wavshap_trl = [];
        
    end
end
