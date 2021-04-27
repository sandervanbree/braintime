function bt_checkallocation(cfg, data, fft_source_constime, fft_source_cutarte)
% This function tests how the two cutmethods (in bt_analyzesources) differ
% in which data they allocate to which cycle. The brain time toolbox
% compresses data cycle-by-cycle after brain time warping, so an important
% quality check is whether cutartefact and consistenttime yield roughly the
% same data-to-cycle allocation. Ideally, differences are relatively small.
% With large differences, it should be noted that both methods could result
% in varying outcomes when used in further analyses.
%
% The function loops through copy of bt_clocktobrain
%
% Use:
% bt_checkallocation(ct_data, fft_source_constime, fft_source_cutarte)
%
% Input Arguments:
% data               % Clock time data structure consisting of both
%                    % classes.
%                    %
%                    % fft_source_constime: toolbox structure obtained from
%                    % bt_selectsource. This data must have
%                    % 'consistenttime' as its cutartefact method (used in
%                    % bt_analyzesources).
%                    %
%                    % fft_source_cutarte: toolbox structure obtained from
%                    % bt_selectsource. This data must have
%                    % 'cutartefact' as its cutartefact method (used in
%                    % bt_analyzesources).
%                    %
% Output:            %
%                    % Figure that demonstrates which cycles
%                    %

%% Obtain warp and phase method
warpmethod = bt_defaultval(config,'warpmethod','stationary');  % Set method for warping (default: stationary)
phasemethod = bt_defaultval(config,'phasemethod','FFT');       % Set phase estimation method used for warping (default: FFT)

%% Check whether input format is correct
if strcmp(fft_source_constime{6},'consistenttime') == 0
    error('The second input is not a structure that used consistenttime during bt_analyzesources');
end

if strcmp(fft_source_cutarte{6},'cutartefact') == 0
    error('The second input is not a structure that used cutartefact during bt_analyzesources');
end

% Transform data structure to raw
if isstruct(data.trial) == 0
    data = ft_checkdata(data,'datatype','raw');
    warning('Timelocked data structure detected, so it was converted to raw');
end

% Set sampling rate of brain time to clock time sampling rate
time   = data.time{1};
phs_sr = round(1/(time(2)-time(1)));

%% Loop over both methods
for method = 1:2
    if method==1
        bt_source = fft_source_constime;
    else
        bt_source = fft_source_cutarte;
    end
    
    % Reset clock time data from any change of previous iteration
    data_temp = data;
    
    % Fetch bt_source information for this particular method
    src_oi = bt_source{1};                                 % Warping source which contains the warping signal
    FFT_phs = cell2mat(bt_source{2});                      % Phase of all frequencies in this warping source
    GED_phs = bt_source{9};                                % Phase of warping signal as estimated using GED
    warpsources = bt_source{3};                            % Warping source data
    srcrank = bt_source{4};                                % Time freq data of selected warping signal
    mintime_fft = bt_source{5}.time(1);                    % Start time of interest
    maxtime_fft = bt_source{5}.time(end);                  % End time of interest
    sr = bt_source{5}.time(2)-bt_source{5}.time(1);        % Sampling rate
    cutmethod = bt_source{6};                              % Applied cutting method
    warpfreq = srcrank(2);                                 % Warped frequency (frequency of the warping signal)
    wvshape = bt_source{7};                                % Average waveshape of the data
    
    
    if strcmp(cutmethod,'cutartefact')                     % Depending on cutmethod, specify original time window of interest
        mintime = mintime_fft+0.5;
        maxtime = maxtime_fft-0.5;
    else
        mintime = mintime_fft;
        maxtime = maxtime_fft;
    end
    
    mintime_ind = nearest(bt_source{5}.time,mintime);    % Index of start time of interest (differs for cutartefact)
    maxtime_ind = nearest(bt_source{5}.time,maxtime);    % Index of end time of interest
    
    %% Warp component's data to template phase vector (based on power oscillation)
    bt_data     = data_temp;
    nsec        = bt_data.time{1}(end)-bt_data.time{1}(1);         % number of seconds in the data
    Ncycles     = warpfreq*nsec;                                   % number of cycles
    cycledur    = round(phs_sr*nsec/Ncycles);                      % samples for cycle
    tempsr      = Ncycles*cycledur/nsec;
    timephs     = linspace(0,Ncycles,phs_sr*nsec);                 % time vector of the unwrapper phase
    ct_phs      = linspace(-pi,(2*pi*Ncycles)-pi,tempsr*nsec);     % set up phase bins for unwrapped phase (angular frequency)
    
    %% Remove the component from original data if desired (default = yes)
    cfg           = [];
    cfg.component = src_oi;
    if isfield(config,'removecomp')
        if strcmp(config.removecomp,'yes')
            if isfield(warpsources,'cfg') && isfield(warpsources.cfg,'method')
                if strcmp(warpsources.cfg.method,'runica')
                    data_temp = ft_rejectcomponent (cfg, warpsources, data_temp);
                else
                    error(['It is not possible to remove the warping component, because'...
                        ' the warping sources are not ICA components.']);
                end
            else
                error(['It is not possible to remove the warping component, because'...
                    ' the warping sources are not ICA components.']);
            end
        end
    else
        % if the removal option is not specified, remove by default
        if isfield(warpsources,'cfg') && isfield(warpsources.cfg,'method') && strcmp(warpsources.cfg.method,'runica')
            data_temp = ft_rejectcomponent (cfg, warpsources, data_temp);
        end
    end
    
    %% Cut out the time window of phase estimation used during bt_analyzesources
    cfg        = [];
    if strcmp(cutmethod,'cutartefact')
        cyclesample = round((1/warpfreq)*1/sr); % Calculate how many samples one cycle consists of
        cfg.toilim = [mintime_fft+0.5-(1/warpfreq) maxtime_fft-0.5+(1/warpfreq)]; % Cut to the time window of interest, plus one cycle
        FFT_phs = FFT_phs(:,mintime_ind-cyclesample:maxtime_ind+cyclesample); % Cut to the time window of interest, plus one cycle
        GED_phs = GED_phs(:,mintime_ind-cyclesample:maxtime_ind+cyclesample); % Do the same for GED phase
    elseif strcmp(cutmethod,'consistenttime')
        cfg.toilim = [mintime_fft maxtime_fft];
    end
    data_temp       = ft_redefinetrial(cfg, data_temp);
    
    %% Opt for FFT or GED phase and adapt data
    if strcmp(phasemethod,'FFT')
        phs = FFT_phs;
    elseif strcmp(phasemethod,'GED')
        
        % Sanity check whether the two phase vectors are within 5% of each other's length
        if abs(size(FFT_phs,2)-size(GED_phs,2)) > 0.05*max(size(FFT_phs,2),size(GED_phs,2))
            error(['The length of the GED estimated phase substantially differs from '...
                'the length of the FFT estimated phase. Please change to config.phasemethod '...
                '= ''FFT'' or change the time window tested during bt_analyzechannels.'])
        else
            phs = GED_phs;
        end
    end
    
    % Check whether phase and data are of the same length
    if abs(length(time)-size(phs,2))>1
        warning(['The phase and data length in length by more than 1 samples. This'...
            ' may mean the wrong parameters were used during bt_analyzesources.']);
    elseif abs(length(time)-size(phs,2))>10
        error(['The phase and data length in length by more than 10 samples. This'...
            ' may mean the wrong parameters were used during bt_analyzesources.']);
    end
    
    % Do some additional slicing to fix slight differences
    if size(phs,2) > length(time)
        phs=phs(:,1:length(time));
    elseif length(time) > size(phs,2)
        cfg       = [];
        cfg.toilim = [time(1) time(end)];
        data_temp      = ft_redefinetrial(cfg, data_temp);
    end
    
    %% Warp the phase of the warping signal to the phase of a stationary oscillation
    bt_data     = data_temp;
    nsec        = bt_data.time{1}(end)-bt_data.time{1}(1);        % number of seconds in the data
    Ncycles     = warpfreq*nsec;                                  % number of cycles
    cycledur    = round(phs_sr*nsec/Ncycles);                     % samples for cycle
    tempsr      = Ncycles*cycledur/nsec;
    timephs     = linspace(0,Ncycles,phs_sr*nsec);                % time vector of the unwrapper phase
    
    if strcmp(warpmethod,'stationary')                            % warp using stationary sinusoid
        ct_phs = linspace(-pi,(2*pi*Ncycles)-pi,tempsr*nsec);     % set up phase bins for unwrapped phase (angular frequency)
        
    elseif strcmp(warpmethod,'waveshape')                         % warp using average waveshape
        wvshape_sm = smoothdata(wvshape,'gaussian',25);           % smooth substantially - 25 bins seems the sweetspot for 201 bins
        [~,trgh] = findpeaks(-wvshape_sm);                        % find troughs
        
        if numel(trgh)~=2                                         % if there are not 2 troughs, this likely
            error(['The waveshape of the warping signal is',...   % means the data is too noisy
                ' too noisy. Please select cfg.method =',...
                ' ''stationary'' and try again.']);
        end
        
        waveshape_cut = wvshape_sm(trgh(1)+1:trgh(2));            % cut waveshape to one cycle (-cos)
        tempsig = repmat(waveshape_cut,[1 floor(Ncycles)]);       % repeat wave Ncycles times (round down)
        misscycs = Ncycles-floor(Ncycles);                        % check if rounding down lost us anything
        extrasamps = round(numel(waveshape_cut)*misscycs);        % how many samples of the cut cycle were lost
        
        if misscycs~=0                                            % append those to the template signal
            tempsig = [tempsig, waveshape_cut(1:extrasamps)];
        end
        
        tempsig_hb = angle(hilbert(tempsig));                     % get the phase using Hilbert transform
        ct_phs = unwrap(tempsig_hb);                              % unwrap the phase
        ct_phs = imresize(ct_phs,[1 tempsr*nsec]);                % resize to the desired length
    end
    
    %% Brain time warp both methods, and save their output
    for nt=1:size(phs,1)
        bt_phs = unwrap(phs(nt,:));
        % Warp phase of single trial onto template phase
        [~,ix,iy] = dtw(bt_phs,ct_phs);
        
        % Find the start of each cycle in the clock to brain time path
        cycles        = false(1,length(iy));
        cycles(1,end) = true;
        for p = 1:length(iy)-1            % Loop through path
            if rem(iy(p),cycledur) == 0   % If we have reached cycledur
                if iy(p) ~= iy(p+1)       % And the next bin is not a repetition
                    cycles(p) = true;     % Label the start of the cycle
                end
            end
        end
        [~, c]=find(cycles == true);      % Find the indices of cycle starts
        xtime=data_temp.time{1}(ix);      % Warped clock time vector
        
        
        % Resample the first cycle's data based on the warping path from brain
        % to clock time
        % First cycle:
        cyl   = bt_data.trial{1,nt}(:,ix(1:c(1)));
        tmpcy = imresize(cyl,[size(bt_data.label,1) cycledur]);  % Resize to cycle duration
        tmptrl(:,1:cycledur) = tmpcy;                            % Replace data
        
        %Remaining cycles
        for cy    = 2:round(Ncycles)
            cyl   = bt_data.trial{1,nt}(:,ix(c(cy-1)+1:c(cy)));
            tmpcy = imresize(cyl,[size(bt_data.label,1) cycledur]);
            tmptrl(:,(cy-1)*cycledur+1:cy*cycledur) = tmpcy;
        end
        
        bt_data.time{1,nt}=timephs;
        
        iycyc=zeros(1,length(ix)); %iycyc = cycle vector for warped phase, using the indices of c
        iycyc(1:c(1,1))=1;
        for cy=2:length(c)
            iycyc(c(cy-1)+1:c(cy))=cy;
        end
        
        % 1st var = indices of the warping trial phase
        % 2nd var = clock time point corresponding to each value of ix
        % 3rd var = indices of the warping template phase
        % 4th var = cycles vector
        x=[ix';xtime;iy';iycyc];
        
        if method == 1 % Save warping information
            consistenttime{1,nt} = x;
        else
            cutartefact{1,nt} = x;
        end
    end
end

%% Check data-to-cycle allocation between the two methods
figure; hold on;

ntr = length(cutartefact); %number of trials
ncyc = cutartefact{1}(4,end); % number of cycles
extracyc=1; % We have added one extra cycle at the beginning and the end to avoid the artefact in the first cycle (from 0s). So, we are cutting out this number of cycles later (this time points are not in the consistent time method)
timepoints=consistenttime{1}(2,:);
timepoints=unique(timepoints); %vector from 0 to 1 second without repetitions
coinc=ones(ntr,length(timepoints)); %coinc= time point is classified in the same cycle or not
coinc=coinc/2; % 0.5 = NaN (to visualize the nans in the final plot)
for tr=1:ntr % trials loop
    
    cycle=cutartefact{tr}; %cycle= cutartefact for the current trial
    % this loop change the 4th row (cycles) from whole numbers to decimal numbers
    for cy=1:cycle(4,end)
        [~, z]=find(cycle(4,:)==cy);
        cycvector=linspace(cy,cy+1,length(z)+1);
        cycle(4,z)=cycvector(1:end-1);
    end
    
    % get the index to cut out the extra cycles (added only in the cut
    % method) and reduce to the warped frequency cycles
    if rem(extracyc,1)==0
        [stidx]=nearest(cycle(4,:),2);
        [endidx]=nearest(cycle(4,:),Ncycles+1);
    end
    
    %if we got several indixes in the previous loop, take only one
    if length(stidx)==2
        stidx=stidx(1,2);
    elseif  length(endidx)==2
        endidx=endidx(1,2);
    end
    cycle=cycle(:,stidx:endidx); % cut the data
    
    % grab the clock time point in which each trial starts/ends
    longitude(tr,1)=cycle(2,1); %start
    longitude(tr,2)=cycle(2,end); %end
    %with this variable we can check the start and end of the trials
    %(and not only its length)
    
    %rename the cycle, from 1 to the warpping frequency (& no more decimal
    %numbers)
    a=round(cycle(4,1),1);
    for tp=1:length(cycle)
        for cy=a:1:cycle(4,end)
            if cycle(4,tp)>=cy && cycle(4,tp)<=cy+1
                cycle(4,tp)=cy-(a-1);
            end
        end
    end
    
    % check if the same clock time point is in the same cycle in long
    % VS short data
    %1st time point
    [~, c_cy]=find(cycle(2,:)>0); %find the first value bigger than 0 in long data (since in short data, the 1st time point is 0)
    [c]= nearest (timepoints,cycle(2,c_cy(1,1))); %use the previous index to find the time in the clock vector time
    [c_short]=nearest(consistenttime{tr}(2,:),cycle(2,c_cy(1,1))); %find when short data is the same as long data
    
    %grab the equivalence (or difference) in a new variable --> coinc. This
    %variable length is the same as clock time vector (as we want to know
    %in each time point is assigned to the same cycle for each method)
    coinc(tr,c)=cycle(4,c_cy(1,1))==consistenttime{tr}(4,c_short(1,1)); %is the same clock time point in the same cycle? Yes --> 1; no --> 0
    
    %Remaining points
    ct=c+1; %counter to put in the correct index of coinc the coincidence value
    for tp=c_cy(1,1)+1:length(cycle)
        if cycle(2,tp)>timepoints(end) % if the last time point of long data is smaller than 1 second, let the rest of the coincidence values be 0.5
        elseif cycle(2,tp)~=cycle(2,tp-1) % if it is the same clock time, don't do it again
            c = nearest(consistenttime{tr}(2,:),cycle(2,tp)); % find when short data is the same as long data
            coinc(tr,ct)=cycle(4,tp) == consistenttime{tr}(4,c(1,1));
            ct=ct+1;
        end
    end
    
end

% here, we can see if each clock time point is assigned to the same cycle
% by both methods. I guess it is not a problem if the same time points
% group is assigned to the same cycle, but it would be if different parts
% if the group were assigned to different cycles. Is it worth changing the
% script to achieve that?

imagesc(coinc)
ylabel('Trials')
xlabel('Time')
% plot interpretation:
% -yellow    = both methods assign the same cycle;
% -blue      = each method assigns a different cycle;
% -turquoise = cut artifact doesnt use these data (warping starts before 0 sec)

% plot start/end trials
figure,

subplot(3,7,1)
plot(longitude{1},'.')
ylabel('Clock Time')
xlabel('Trials')
ylim([-.5 1.5])