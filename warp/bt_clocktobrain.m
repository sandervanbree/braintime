function [bt_struc] = bt_clocktobrain(config, data, bt_comp)
% Help to be added

compoi = bt_comp{1};
phs = bt_comp{2};
comp = bt_comp{3};
topcomps = bt_comp{4};
mintime = bt_comp{5}.time(1);
maxtime = bt_comp{5}.time(end);
cutmethod = bt_comp{6};
warpfreq = topcomps(2);

if isfield(config,'btsrate')
    phs_sr = config.btsrate;
else
    phs_sr = 512; %Default sampling rate
end

%remove the component from original data
cfg           = [];
cfg.component = compoi;
data          = ft_rejectcomponent (cfg, comp, data);

% cut out the time window of fft (from which the phase was extracted)
cfg        = [];
cfg.toilim = [mintime maxtime];
data       = ft_redefinetrial(cfg, data);

%% Step 8: Warp component's data to template phase vector (based on power oscillation)
% Re-Organize EEG data by phase
bt_data=data;
nsec=maxtime-mintime;
Ncycles_pre=warpfreq*nsec; %number of cycles * seconds
cycledur=round(phs_sr*nsec/Ncycles_pre); %samples for cycle
tmp_sr=Ncycles_pre*cycledur/nsec;
tempphs=linspace(-pi,(2*pi*Ncycles_pre)-pi,tmp_sr*nsec);% set up phase bins for unwrapped phase (angular frequency)
timephs=linspace(0,Ncycles_pre,phs_sr*nsec); %time vector of the unwrapper phase

for nt=1:size(phs{1},1)
    tmpphstrl=unwrap(phs{1}(nt,:));
    % Warp phase of single trial onto template phase
    [~,ix,iy] = dtw(tmpphstrl,tempphs); %to get the equivalence index between template and trial phase
    
    %how long is each cycle?
    cycles=zeros(1,length(iy));
    for tp=2:length(iy)
        if rem(iy(tp),cycledur)==0
            if iy(tp)~=iy(tp-1)
                cycles(tp)=500; %label each cycle start with arbitrary number
            end
        end
    end
    [~, c]=find(cycles==500); 
    
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
    
    % Create warped trials
    bt_data.trial{1,nt}=imresize(tmptrl,[size(tmptrl,1) phs_sr*nsec]);
    bt_data.time{1,nt}=timephs;
end

% If method is cut artefact, cut right time window and adjust time vector
if strcmp(cutmethod,'cutartefact')
    % correct time window of interest
    mintime = mintime+0.5;
    maxtime = maxtime-0.5;
    
    % cut out the extra cycles
    FreqDiff = Ncycles_pre-(warpfreq*0.5)-warpfreq*((maxtime)-(mintime)); %How much is the time vector off?
    
    startind = (phs_sr*0.5)+1; % first cycle of the time window of interest index
    endind   = length(bt_data.time{1,1})-(phs_sr*0.5); %What's the index of the end of the desired new data?
    
    cfg         = [];
    cfg.latency = [bt_data.time{1}(startind) bt_data.time{1}(endind)];
    bt_data  = ft_selectdata(cfg,bt_data);
    bt_data.trialinfo = data.trialinfo;
    
    % correct the cycles vector
    for trl = 1:numel(bt_data.trial)
        bt_data.time{trl} = bt_data.time{trl}-FreqDiff;
    end
end

% reformat data structure and include basic information
bt_struc.data = bt_data;
bt_struc.toi = [mintime maxtime];
bt_struc.freq = warpfreq; %Warped frequency
bt_struc.clabel = bt_data.trialinfo;
