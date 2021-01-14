function [bt_struc] = bt_clocktobrain(config, data, bt_carrier)
% #### SCRIPT WORKING, LINGO UPDATED SOON ####
% Warp clock to brain time. The clock time data is resampled based on the
% warping path from the brain time phase vector to the phase of a
% stationary sinusoid.
%
% Use:
% [bt_struc] = bt_clocktobrain(config,data,bt_carrier)
%
% Input Arguments:
% config
%   - btsrate        % Sampling rate of the brain time data.
%                    %
%   - removecomp     % 'yes': removes component from the brain time data.
%                    % When analyzing brain time data using your own
%                    % analysis pipeline, you may wish to remove the 
%                    % component to avoid circularity. See the brain time
%                    % toolbox paper for more details.
%                    %
%                    % 'no': keeps component in the brain time data.
%                    %
% data               % Preprocessed clock time data structure consisting of
%                    % both classes.
%                    %
% bt_carrier         % Data structure obtained from bt_choosecarrier.
%                    % Includes: Carrier's time frequency information, and 
%                    % config details saved for later retrieval.
%                    %
% Output:            %
% bt_struc           % Data structure with: brain time data, its 
%                    % time frequency information, and config details
%                    % saved for later retrieval.

%% Get basic info
channeloi = bt_carrier{1};                            % Channel which contains the carrier
phs = cell2mat(bt_carrier{2});                        % Phase of all frequencies in this channel
channels = bt_carrier{3};                             % Channel data 
chanrank = bt_carrier{4};                             % Time freq data of chosen channel
mintime = bt_carrier{5}.time(1);                      % Start time of interest
maxtime = bt_carrier{5}.time(end);                    % End time of interest
sr = bt_carrier{5}.time(2)-bt_carrier{5}.time(1);     % Sampling rate
cutmethod = bt_carrier{6};                            % Applied cutting method 
warpfreq = chanrank(2);                               % Warped frequency (frequency of the carrier)
mintime_ind = bt_carrier{7}(1);                       % Index of start time of interest (differs for cutartefact) 
maxtime_ind = bt_carrier{7}(2);                       % Index of end time of interest (differs for cutartefact)
analyzemethod = bt_carrier{8};                        % Save carrier choosing method 

% Set up sampling rate
if isfield(config,'btsrate')
    phs_sr = config.btsrate;
else
    phs_sr = 512; %Default sampling rate
end

%% Remove the component from original data if desired (default = yes)
cfg           = [];
cfg.component = channeloi;
if isfield(config,'removecomp')
    if strcmp(config.removecomp,'yes')
        if strcmp(analyzemethod,'bt_choosecarrier')
        data = ft_rejectcomponent (cfg, channels, data);
        else
        error('For component analysis using bt_GEDanalyzechoose, you can remove the component by specifying cfg.removecomp to ''yes'' in bt_GEDanalyzechoose.')
        end
    end
else
    % if the removal option is not specified, remove by default
    if isfield(channels.cfg,'method')&&strcmp(analyzemethod,'bt_choosecarrier') %only makes sense if channel data are ICA components and if method was bt_chooseanalyze
        if strcmp(channels.cfg.method,'runica')
    data = ft_rejectcomponent (cfg, channels, data);
        end
    end
end

%% Cut out the time window of fft (from which the phase was extracted)
cfg        = [];
if strcmp(cutmethod,'cutartefact')
    cyclesample = round((1/warpfreq)*1/sr); %Calculate how many samples one cycle consists of
    cfg.toilim = [mintime+0.5-(1/warpfreq) maxtime-0.5+(1/warpfreq)]; %Cut to the time window of interest, plus one cycle
    phs = phs(:,mintime_ind-cyclesample:maxtime_ind+cyclesample); %Cut to the time window of interest, plus one cycle
elseif strcmp(cutmethod,'consistenttime')
    cfg.toilim = [mintime maxtime];
end
data       = ft_redefinetrial(cfg, data);

% Check whether phase and data are of the same length
    if abs(length(data.time{1})-length(phs))>1
        warning('Warning: phase vector and data differ in length by more than 1 sample.');
    end   
if length(phs) > length(data.time{1})
    phs=phs(:,1:length(data.time{1}));
elseif length(data.time{1}) > length(phs)
    cfg       = [];
    cfg.toilim = [data.time{1}(1,1) data.time{1}(1,length(phs))];
    data      = ft_redefinetrial(cfg, data);
end

%% Warp component's data to template phase vector (based on power oscillation)
% Re-Organize EEG data by phase
bt_data=data;
nsec=bt_data.time{1}(end)-bt_data.time{1}(1); %number of seconds in the data
Ncycles_pre=warpfreq*nsec; %number of cycles * seconds
cycledur=round(phs_sr*nsec/Ncycles_pre); %samples for cycle
tmp_sr=Ncycles_pre*cycledur/nsec;
tempphs=linspace(-pi,(2*pi*Ncycles_pre)-pi,tmp_sr*nsec);% set up phase bins for unwrapped phase (angular frequency)
timephs=linspace(0,Ncycles_pre,phs_sr*nsec); %time vector of the unwrapper phase

for nt=1:size(phs,1)
    tmpphstrl=unwrap(phs(nt,:));
    % Warp phase of single trial onto template phase
    [~,ix,iy] = dtw(tmpphstrl,tempphs); %to get the equivalence index between template and trial phase
    
    %how long is each cycle?
    cycles=zeros(1,length(iy));
    cycles(1,end)=500;
    for tp=1:length(iy)-1
        if rem(iy(tp),cycledur)==0
            if iy(tp)~=iy(tp+1)
                cycles(tp)=500; %label the start of each cycle with '500'
            end
        end
    end
    [~, c]=find(cycles==500); %find the indices of the starts
    
    %get equal samples by cycle
    %First cycle
    cyl=bt_data.trial{1,nt}(:,ix(1:c(1)));
    tmpcy=imresize(cyl,[size(bt_data.label,1) cycledur]);
    tmptrl(:,1:cycledur)=tmpcy;
    
    %Remaining cycles
    for cy=2:Ncycles_pre
        cyl=bt_data.trial{1,nt}(:,ix(c(cy-1)+1:c(cy)));
        tmpcy=imresize(cyl,[size(bt_data.label,1) cycledur]);
        tmptrl(:,(cy-1)*cycledur+1:cy*cycledur)=tmpcy;
    end
    
    %Create warped trials
    bt_data.trial{1,nt}=imresize(tmptrl,[size(tmptrl,1) numel(timephs)]);
    bt_data.time{1,nt}=timephs;
end

% If method is cut artefact, cut right time window and adjust time vector
if strcmp(cutmethod,'cutartefact')
    % correct time window of interest to true time
    mintime = mintime+0.5;
    maxtime = maxtime-0.5;
    
    startind = nearest(bt_data.time{1},1); %Find index of first cycle in window of interest
    endind = nearest(bt_data.time{1},warpfreq*(maxtime-mintime)+1); %Find index of last cycle in window of interest
    
    cfg         = [];
    cfg.latency = [bt_data.time{1}(startind) bt_data.time{1}(endind)]; % Cut to time window of interest
    bt_data  = ft_selectdata(cfg,bt_data);
    bt_data.trialinfo = data.trialinfo;
    
    % correct the cycles vector
    for trl = 1:numel(bt_data.trial)
        bt_data.time{trl} = bt_data.time{trl}-1; %cycles vector is off by 1
    end
end

%% Reformat data structure and include basic info
bt_struc.data = bt_data;                                     % Brain time warped data
bt_struc.toi = [mintime maxtime];                            % Start and end time of interest
bt_struc.freq = warpfreq;                                    % Warped frequency (frequency of the carrier)
bt_struc.clabel = bt_data.trialinfo;                         % Classification labels
