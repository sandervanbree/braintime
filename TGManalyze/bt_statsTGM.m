function bt_statsTGM(config, bt_data, bt_quantTGM)
% Test whether quantified recurrence in the TGM is statistically reliable.
% Performs one or two-level permutation testing, where null distributions
% are created by shuffling the classification labels nperm1(*N) times and
% collecting the power spectra from the resulting AC maps. The power
% spectra of the empirical data are compared against this null
% distribution.
%
% Use:
% bt_statsTGM(config, bt_data, bt_quantTGM)

% Establish basic parameters
toi = bt_quantTGM.toi;
warpfreq = bt_quantTGM.warpfreq;
acfft = bt_quantTGM.acfft;
modefreq = bt_quantTGM.modefreq;
modefreqind = bt_quantTGM.modefreqind;
TGM = bt_quantTGM.TGM;
refdimension = bt_quantTGM.refdimension;
timevec = bt_quantTGM.timevec;
normalizer = bt_quantTGM.normalizer;
clabel = config.clabel;
cfg_mv = config.mvpacfg;

if isfield(config,'statsrange') %What is the recurrence rate range of interest?
    statsrange = config.statsrange(1):config.statsrange(2);
else
    statsrange = 1:40;
end

% Calculate autocorrelation map (AC)
ac=autocorr2d(TGM);

% Run FFT over all rows and columns of the AC map
nvecs=numel(ac(:,1));

for vec=1:nvecs
    %1st dimenssion
    [PS,f]=Powspek(ac(vec,:),nvecs/normalizer);
    PS1(vec,:) = PS(statsrange);
    
    %2nd dimension
    [PS,f]=Powspek(ac(:,vec),nvecs/normalizer);
    PS2(vec,:) = PS(statsrange);
    
end
avg_PS = mean(PS1,1)+mean(PS2,1); %Mean power spectra
modepow_emp = avg_PS(modefreqind); %Mean power at mode freq
fullspec_emp = avg_PS;

if config.permlevels == 1
    for perm = 1:config.numperms1
        fprintf('First level permutation number %i\n', perm);
        clabel = clabel(randperm(numel(clabel)));
        [permTGM, ~] = mv_classify_timextime(cfg_mv, bt_data.trial, clabel);
        
        % Calculate autocorrelation map (AC)
        ac=autocorr2d(permTGM);
        
        % Run FFT over all rows and columns of the AC map
        nvecs=numel(ac(:,1));
        
        for vec=1:nvecs
            %1st dimenssion
            [PS,f]=Powspek(ac(vec,:),nvecs/normalizer);
            PS1(vec,:) = PS(statsrange);
            
            %2nd dimension
            [PS,f]=Powspek(ac(:,vec),nvecs/normalizer);
            PS2(vec,:) = PS(statsrange);
            
        end
        avg_PS = mean(PS1,1)+mean(PS2,1); %Mean power spectra
        modepow_shuff(perm) = avg_PS(modefreqind); %Mean power spectra at mode freq
        fullspec_shuff(perm,:) = avg_PS;
    end
    mean_modepow_shuff = mean(modepow_shuff);
end

f=f(statsrange); %filter frequency vector based on range of interest

%create confidence interval for each frequency bin
for freq = 1:numel(f)
    
    %Grab data
    temp_data = fullspec_shuff(:,freq);
    
    %Calculate mean
    temp_mean = mean(temp_data);
    
    %Calculate standard mean error (SEM)
    temp_SEM = std(temp_data)/sqrt(length(temp_data));
    
    %Calculate t-score
    temp_ts = tinv([0.025  0.975],length(temp_data)-1);
    
    %Calculate confidence interval (CI)
    temp_CI = temp_mean + temp_ts*temp_SEM;
    
    freq_CI(freq,:) = temp_CI; % save CI
end

% Plot results
% 1st plot: relationship empirical and shuffled amp at mode freq
figure
subplot(1,2,1)
x = ones(1,numel(modepow_shuff));
plot(x,modepow_shuff,'k*','MarkerSize',5,'LineWidth',3);
grid on;
hold on;
plot(x(1),modepow_emp,'r*','MarkerSize',15,'LineWidth',3);
legend('Shuffled mode power','Empirical mode power')
ax = gca;
ax.XTick = 1;
ax.XTickLabels = {[' ']};
% Set up axes.
ylabel('Mean amplitude at mode frequency')
xlabel(sprintf('%.2f Hz',modefreq))
xlim([0.8, 1.2]);
title(sprintf('Amplitude of mode frequency (%.2f Hz)',modefreq))

% 2nd plot: relationship empirical and shuffled amp at all freq
subplot(1,2,2); hold on
p1 = plot(f,fullspec_emp,'LineWidth',2,'Color','b');
hold on
low_CI = freq_CI(:,1)';
hi_CI = freq_CI(:,2)';
p2 = plot(f,low_CI,'LineWidth',2,'Color','k');
p3 = plot(f,hi_CI,'LineWidth',2,'Color','k');

if strcmp(refdimension,'braintime') %warp freq line is dependent on clock (warped freq) or brain time (1 hz)
p4 = line([1 1], [0 max(fullspec_emp)],'color',[1 0 1],'LineWidth',2); %Line at warped freq
else
p4 = line([warpfreq warpfreq], [0 max(fullspec_emp)],'color',[1 0 1],'LineWidth',2); %Line at warped freq
end
% set up axes
xlabel('Recurrence frequency')
ylabel('Mean amplitude')
title(sprintf('Amplitude across recurrence rates'))
% cipatch = [low_CI', hi_CI'];
% area(cipatch)
cipatch = fill([f,f], [low_CI, hi_CI], [220 220 220]./255);
set(cipatch, 'edgecolor', 'none');
set(cipatch, 'FaceAlpha', 0.6);
legend([p1 p2 p3 p4],{'Empirical amplitude','Lower limit CI','Higher limit CI','Warped frequency'});


