function bt_checkallocation(config, data, bt_carrier_cons, bt_carrier_cut)

% Set up sampling rate
if isfield(config,'btsrate')
    phs_sr = config.btsrate;
else
    phs_sr = 512; %Default sampling rate
end

for method = 1:2
    if method==1
        bt_carrier = bt_carrier_cons;
    else
        bt_carrier = bt_carrier_cut;
    end
    
    channeloi = bt_carrier{1}; %channel of interest
    phs = cell2mat(bt_carrier{2}); %its phase
    channels = bt_carrier{3}; %channel structure from FieldTrip
    topchans = bt_carrier{4}; %top components
    mintime = bt_carrier{5}.time(1);
    maxtime = bt_carrier{5}.time(end);
    cutmethod = bt_carrier{6};
    warpfreq = topchans(2); %warped frequency
    mintime_ind = bt_carrier{7}(1);
    maxtime_ind = bt_carrier{7}(2);
    
    %% Cut out the time window of fft (from which the phase was extracted)
    cfg        = [];
    if method == 2
        cyclesample = round((1/warpfreq)*channels.fsample); %Calculate how many samples one cycle consists of
        cfg.toilim = [mintime+0.5-(1/warpfreq) maxtime-0.5+(1/warpfreq)]; %Cut to the time window of interest, plus one cycle
        phs = phs(:,mintime_ind-cyclesample:maxtime_ind+cyclesample); %Cut to the time window of interest, plus one cycle
    else
        cfg.toilim = [mintime maxtime];
    end
    data_temp       = ft_redefinetrial(cfg, data);
    
    %% Warp component's data to template phase vector (based on power oscillation)
    % Re-Organize EEG data by phase
    bt_data=data_temp;
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
                    cycles(tp)=500;
                end
            end
        end
        [~, c]=find(cycles==500);
        tx=data.time{1}(ix); % what time point corresponds to each point of the warping trial phase? tx= time vector for warped trial phase
        
        
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
        
       
        bt_data.time{1,nt}=timephs;
        
        iycyc=zeros(1,length(ix)); %iycyc= cycle vector for warped phase, using the indices of c
        iycyc(1:c(1,1))=1;
        for cy=2:length(c)
            iycyc(c(cy-1)+1:c(cy))=cy;
        end
        
        x=[ix';tx;iy';iycyc]; % 1st row= indices of the warping trial phase. 2nd row= clock time point corresponding to each value of ix. 3th row= indices of the warping template phase. 4th row= cycles vector
        
        if method==1
            consistenttime{1,nt} = x;
        else
            cutartefact{1,nt} = x;
            bt_data.time{1,nt}=timephs;
            startind = findnearest(bt_data.time{1},1); %Find index of first cycle in window of interest
            endind = findnearest(bt_data.time{1},warpfreq*(maxtime-mintime)+1); %Find index of last cycle in window of interest
        end
    end
end

%% Check data-to-cycle allocation between the two methods
figure

ntr=length(cutartefact); %number of trials

ncyc=cutartefact{1}(4,end); % number of cycles
extracyc=1; % extra cycles (this number of cycles has to be cut out at the beggining and at the end of the trials)
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
    
    % get the index to cut out the extra cycles
    if rem(extracyc,1)==0
        [~, stidx]=findnearest(cycle(4,:),2);
        [~, endidx]=findnearest(cycle(4,:),warpfreq+1);
    end
    
    if length(stidx)==2
        stidx=stidx(1,2);
    elseif  length(endidx)==2
        endidx=endidx(1,2);
    end
    cycle=cycle(:,stidx:endidx);
    
    % create a variable to know the clock time point in which each
    % trial:
    longitud(tr,1)=cycle(2,1); %starts
    longitud(tr,2)=cycle(2,end); %ends
    %with this variable we can check the start and end of the trials
    %(and not only its length)
    
    
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
    [~, c_cy]=find(cycle(2,:)>0); %find the first value bigger than 0 in long data
    [~, c]=find(timepoints==cycle(2,c_cy(1,1))); %use the previous index to know in which place of time point it should start (we are not using the values before 0)
    [~, c_short]=find(consistenttime{tr}(2,:)==cycle(2,c_cy(1,1))); %find when short data is in the same clock time point as long data
    
    coinc(tr,c)=cycle(4,c_cy(1,1))==consistenttime{tr}(4,c_short(1,1)); %is the same clock time point in the same cycle?
    %loop for the rest of the points
    ct=c+1; %counter to put in the correct index of coinc the coincidence value
    for tp=c_cy(1,1)+1:length(cycle)
        if cycle(2,tp)>timepoints(end) %if the last time point of long data is smaller than 1 second, let the rest of the coincidence values be 0.5 (NaNs)
        elseif cycle(2,tp)~=cycle(2,tp-1) %if it is the same clock time, don´t do it again
            [~, c]=find(consistenttime{tr}==cycle(2,tp));
            coinc(tr,ct)=cycle(4,tp)==consistenttime{tr}(4,c(1,1));
            ct=ct+1;
        end
    end
    
end


subplot(3,7,1),imagesc(coinc)
ylabel('Trials')
xlabel('Time')


% plot start/end trials
figure,

subplot(3,7,1)
plot(longitud{1},'.')
ylabel('Clock Time')
xlabel('Trials')
ylim([-.5 1.5])
