function [asymmidx,asymmidx_t,wavshap] = bt_calcwaveshape(wsources,ncycles,srcrank)
% Written by Simon Hanslmayr. Modified for the Brain Time Toolbox by Sander
% van Bree
% We normalizes the asymmetry index to give a ratio of how much of the
% cylce is covered by ascending vs descending flank using a value which
% ranges between -1 and 1.
% Formula: (t ascending flank - t descending flank)/half a cycle
% Belluscio et al. 2012 J Neurosci

ntrls          = size(wsources.trial,2);                  % Fetch number of trials
totsamples     = numel(wsources.time{1,1});               % Total number of samples per trial
ntopwsources   = numel(srcrank(:,2));                     % Fetch number of top warping sources
wfreqs         = srcrank(:,2);                            % Warping frequency per component
source_inds    = srcrank(:,1);                            % Index (name) of top warping source

% winwidth = 3/centerfreq; % Normalize the filter window to scale with frequency (lower frequencies, larger window)
winwidth = 2;
    

%% Prepare basic variables
% Compute basic parameters
deltat   = 1000/wsources.fsample;
tvec     = -ncycles/2:0.01:ncycles/2;                % Set up time vector for average waveshape

% Pre-allocate primary variables
asc_tmp = [];                                        % Ascending flank of cycles (trough to peak)
desc_tmp = [];                                       % Descending flank of cycles (peak to trough)
wavshap_trl = [];                                    % The average wave shape for a trial
wavshap = zeros(ntopwsources,numel(tvec));           % The average wave shape across all trials
% asymmidx = cell(ntopwsources,1);                     % Asymmetry index from -1 to 1, representing asymmetry direction (ascending or descending flank)

%% Loop over data
    
    % For every warping source
    for w = 1:ntopwsources
        c = 0;                                       % Counter
        
        % Get the warping frequency
        centerfreq = wfreqs(w);
        
        % Band pass filter component data
        cfg = [];
        cfg.channel  = source_inds(w);
        data_curr    = ft_preprocessing(cfg,wsources);      % Extract current warping source (raw)
        
        cfg.bpfilter = 'yes';
        cfg.bpfreq   = [centerfreq-winwidth/2 centerfreq+winwidth/2]; 
        data_bp_curr = ft_preprocessing(cfg,wsources);      % Extract current warping source (band-pass)
        
        % Parameters for this frequency
        win      = round((1000/centerfreq)/deltat)*(ncycles./2); % Window covering 1 cycle
        halfcyc  = ((1000./centerfreq)./2)./deltat;              % Window covering 1/2 cycle
        locwin   = -floor(halfcyc./4):1:ceil(halfcyc./4);        % Window over which local minima and maxima are calculated
        
        % Update user
        disp(['Calculating waveshape and asymmetry metrics for top warping source ',num2str(w),' out of ',num2str(ntopwsources),'...']);
        
        % For every trial
        for n=1:ntrls
            trlf   = data_bp_curr.trial{1,n};          % Extract current trial (band pass filtered)
            trlr   = data_curr.trial{1,n};             % Extract current trial (raw)
            
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

            % Average 
            curr_asymm = (asc_tmp-desc_tmp)./(halfcyc*2);
            asymmidx_trl(n) = mean(curr_asymm);
            
            [~,~,~,tres] = ttest(curr_asymm,0,'Tail','both');
            asymmidx_trl_t(n) = tres.tstat;
            
            asc_tmp = [];
            desc_tmp = [];
        end
        
        asymmidx(w,:) = asymmidx_trl;
        asymmidx_t(w,:) = asymmidx_trl_t;
        
        % This calculates the average waveshape and resizes it
        tmp_mn = mean(wavshap_trl,1);
        wavshap(w,:) = imresize(tmp_mn,[1 numel(tvec)]);
        
        % Reset variables
        asc_tmp = [];
        desc_tmp = [];
        wavshap_trl = [];
        asymmidx_trl = [];

    end
end
