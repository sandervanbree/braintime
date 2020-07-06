function [bt_quantTGM] = bt_quantifyTGM(config, TGM)
%help function to be added

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
    subplot(2,2,2)
    pcolor(timevec,timevec,ac);shading interp;title(['Autocorrelation map'])
    caxis([-max(max(abs(ac))) max(max(abs(ac)))])
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

bt_quantTGM.toi = toi;
bt_quantTGM.warpfreq = warpfreq;
bt_quantTGM.acfft = acfft;
bt_quantTGM.timevec = timevec;
bt_quantTGM.normalizer = normalizer;
bt_quantTGM.modefreq = modefreq;
bt_quantTGM.modefreqind = modefreqind;
bt_quantTGM.TGM = TGM;
bt_quantTGM.refdimension = config.refdimension;
