function [bt_quantTGM] = bt_quantifyTGM(config, TGM)
% Quantify the degree of cross-time recurrence in the time generalization
% matrix (TGM). Creates an autocorrelation map (AC map) of the TGM, and
% applies an FFT over each row and column of the AC map. For each row and
% column, the maximum power frequency will be displayed. The mode frequency
% represents the primary recurrence rate present in the TGM.
%
% Use:
% [bt_quantTGM] = bt_quantifyTGM(config, TGM)
%
% Input Arguments:
% config
%   - bt_struc       % brain time data structure as obtained by 
%                    % bt_clocktobrain.
%                    %
%   - refdimension   % 'clocktime': find TGM recurrence as a function of
%                    % clock time seconds in the brain time data.
%                    % 'braintime': find TGM recurrence as a function of
%                    % cycles of the warped frequency in the brain time
%                     %data.
%                    % 
%   - figure         % 'yes' or 'no': display TGM, AC map, and primary
%                    % frequency for all rows and columns in the AC map.
%                    %
% TGM                % TGM obtained by mv_classify_timextime
%                    %
% Output:            %
% bt_quantTGM        % Data structure with: TGM, AC map, FFT information 
%                    % of the AC map, and config details saved for later
%                    % retrieval.

Warp clock to brain time. The clock time data is resampled based on the
% warping path from the brain time phase vector to the phase of a
% stationary sinusoid.
%
% Use:
% [bt_struc] = bt_analyzecomps(config,data,bt_comp)

% Establish basic parameters
toi = config.bt_struc.toi;
warpfreq = config.bt_struc.freq;
duration = toi(2)-toi(1);
if strcmp(config.refdimension,'braintime')
    normalizer = duration*warpfreq; %normalize by cycles in the data
    timevec = config.bt_struc.data.time{1};
elseif strcmp(config.refdimension,'clocktime')
    normalizer = duration; %normalize by seconds in the data
    timevec = linspace(toi(1),toi(2),numel(config.bt_struc.data.time{1}));
end

% Calculate autocorrelation map (AC)
ac=autocorr2d(TGM);

% Run FFT over all rows and columns of the AC map
nvecs=numel(ac(:,1));

for vec=1:nvecs %Only for half is needed
    %1st dimenssion
    [PS,f]=Powspek(ac(vec,:),nvecs/normalizer);
    [pks,locs]=findpeaks(PS);
    maxpk=find(pks==max(pks));
    
    acfft_dim1(1,vec)=pks(maxpk); %What's the amplitude of the peak?
    acfft_dim1(2,vec)=f(locs(maxpk)); %What's the frequency of the peak?
    
    %2nd dimension
    [PS,f]=Powspek(ac(:,vec),nvecs/normalizer);%timevec(end));
    [pks,locs]=findpeaks(PS);
    maxpk=find(pks==max(pks));
    
    acfft_dim2(1,vec)=pks(maxpk); %What's the amplitude of the peak?
    acfft_dim2(2,vec)=f(locs(maxpk)); %What's the frequency of the peak
end

% put together both dimensions
acfft=[acfft_dim1,acfft_dim2];

% get the overall mode
modefreq=mode(acfft(2,:));
[~,modefreqind]=find(abs(modefreq-f)==min(abs(modefreq-f)));

if isfield(config,'figure')
    if strcmp(config.figure,'yes')
        figopt = 1;
    else
        figopt = 0;
    end
else
    figopt = 1; %Default yes
end

if figopt == 1
    % Plot TGM
    figure;
    subplot(2,2,1)
    cfg_plot= [];
    cfg_plot.x   = timevec;
    cfg_plot.y   = cfg_plot.x;
    mv_plot_2D(cfg_plot, TGM);
    cb = colorbar;
    title(cb,'accuracy')
    xlim([timevec(1) timevec(end)]);
    ylim([timevec(1) timevec(end)]);
    xticks(yticks) % make ticks the same on the two axes
    title(['Time Generalization Matrix'])
    if strcmp(config.refdimension,'braintime')
        xlabel('Test data (cycles)')
        ylabel('Training data (cycles)')
    elseif strcmp(config.refdimension,'clocktime')
        xlabel('Test data (seconds)')
        ylabel('Training data (seconds)')
    end
    
    % Plot AC map
    % Detect appropriate color range by zscoring  
    mn_ac=median(ac(:));
    sd_ac=std(ac(:));
    ac_z=(ac-mn_ac)/sd_ac;
    indx = (abs(ac_z)<5); %filter to include only z-scores under 5
    clim = max(ac(indx)); %take the max number as clim for plotting
    
    subplot(2,2,2)
    pcolor(timevec,timevec,ac(1:numel(timevec),1:numel(timevec)));shading interp;title(['Autocorrelation map'])
    caxis([-clim +clim])
    cb=colorbar;
    title(cb,'corr')
    if strcmp(config.refdimension,'braintime')
        xlabel('Shift by x-cycle')
        ylabel('Shift by y-cycle')
    elseif strcmp(config.refdimension,'clocktime')
        xlabel('Shift by x-sec')
        ylabel('Shift by y-sec')
    end
    hold on
    
    % Plot mode frequency
    subplot(2,2,3)
    plot(acfft(2,:),1:nvecs*2,'.')
    title(['Mode peak frequency: ',num2str(modefreq)])
    xlabel('Peak frequency')
    ylabel('AC Row or column')
end

%% Save basic info
bt_quantTGM.toi = toi;
bt_quantTGM.warpfreq = warpfreq;
bt_quantTGM.acfft = acfft;
bt_quantTGM.timevec = timevec;
bt_quantTGM.normalizer = normalizer;
bt_quantTGM.modefreq = modefreq;
bt_quantTGM.modefreqind = modefreqind;
bt_quantTGM.TGM = TGM;
bt_quantTGM.refdimension = config.refdimension;
