%%% In tutorial 9 we will transform perform more commonplace electro-
%%% physiology analyses to show that brain time warped data is also useable
%%% outside of the brain time toolbox.
%%% To do this, we will take the brain time warped data from tutorial 1,
%%% and compare the event-related potentials (ERPs), power spectra, and
%%% intertrial coherence between clock and brain time.
%%%
%%% Each analysis of this script can be carried out without having the
%%% brain time toolbox installed (though the toolbox is necessary to obtain
%%% brain time warped data). When using warped data outside of the toolbox,
%%% it is important to be coignizant of circularity concerns. For details
%%% on these concerns, please refer to brain time documentation on
%%% https://github.com/sandervanbree/braintime

% Load output
load tutorial1_output
load dipolesim_tutorial dipolesim_params

% Display ground truth model
disp(dipolesim_params.mainfreq); % Simulated brain time frequency
disp(dipolesim_params.dipfreqs); % Simulated 1/f random frequencies. These
                                 % frequencies are also in the data, but
                                 % less prominently. They may show up in
                                 % results.
                                 
%% Analysis 1: Event-related potentials (ERPs)
% If brain time warping overcomes nonstationarities in the data, the
% average signal should look more aligned across trials after warping. To
% test this, let us split up the data into left and right hemifield trials
% and compare the ERPs between clock and brain time.
                                  
% Find trials of each condition for clock and brain time data
cfg         = [];
cfg.trials  = find(ct_data.trialinfo == 1); % LEFT
ct_left     = ft_selectdata(cfg,ct_data);
bt_left     = ft_selectdata(cfg,bt_warpeddata.data);
cfg.trials  = find(ct_data.trialinfo == 2); % RIGHT
ct_right    = ft_selectdata(cfg,ct_data);
bt_right   = ft_selectdata(cfg,bt_warpeddata.data);

% Perform a timelocked analysis (average across trials)
cfg      = [];
ct_left  = ft_timelockanalysis(cfg,ct_left);
bt_left  = ft_timelockanalysis(cfg,bt_left);
ct_right = ft_timelockanalysis(cfg,ct_right);
bt_right = ft_timelockanalysis(cfg,bt_right);

% Plot results
cfg = [];
cfg.layout = 'layout_tutorial.mat';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.channel     = 'A15'; % Electrode in the left hemisphere (e.g., 'A15') or 
                         % right hemisphere (e.g., 'A28').

% Plot clock time
ft_singleplotER(cfg, ct_left, ct_right);
title('Clock time ERP');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Left condition','Right condition');

% Plot brain time
ft_singleplotER(cfg, bt_left, bt_right);
title('Brain time ERP');
xlabel('Time (cycles)'); ylabel('Amplitude');
legend('Left condition','Right condition');

% We see increased oscillatory structure in brain time as compared to clock
% time. The simulated anti-phase relation is not apparent here because all
% the dipoles are simulated rather close to each other -- causing the
% signals to sum spatially.
% Results vary across parameters.

%% Analysis 2: Power spectrum
% If brain time warping overcomes nonstationarities in the data, the
% power spectrum should show an enhanced peak at the brain time warped
% frequency after warping.

% Pwelch parameters
segmentLength = 100; % number of samples per segment
noverlap      = 25; % overlapping samples between segments
fs            = ct_data.fsample; % sampling rate

% Clock time
for tr = 1:numel(ct_data.trial) % Loop over trials
    disp([num2str(round(tr./numel(ct_data.trial),2).*100),'% of Clock time power spectrum']);
    for c = 1:numel(ct_data.label) % and over channels
       ct_trial = ct_data.trial{tr}(c,:); 
       [ct_fft(c,tr,:) , freq] = pwelch(ct_trial,segmentLength,noverlap,fs);
    end
end
ct_ps = squeeze(mean(squeeze(mean(ct_fft,1)),1)); % Average across trials and channels

% Brain time
for tr = 1:numel(bt_warpeddata.data.trial)
   disp([num2str(round(tr./numel(bt_warpeddata.data.trial),2).*100),'% of Brain time power spectrum']);
    for c = 1:numel(bt_warpeddata.data.label)
       bt_trial = bt_warpeddata.data.trial{tr}(c,:); 
       [bt_fft(c,tr,:) , freq] = pwelch(bt_trial,segmentLength,noverlap,fs);
    end
end
bt_ps = squeeze(mean(squeeze(mean(bt_fft,1)),1));

% Plot results
figure; hold on; title('Power spectrum');
plot(ct_ps,'LineWidth',3);hold on;
plot(bt_ps,'LineWidth',3);
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([2 30]);
legend({'clock time','brain time'});

% Depending on simulation paramers and the chosen warping signal, we may
% see a stronger peak in brain time than clock time.
% Results vary across parameters.

%% Analysis 3: Intertrial coherence (ITC)
% If brain time warping overcomes nonstationarities in the data, the
% intertrial coherence around the warping frequency should be enhanced
% after warping.

toi    = 1:0.02:2;        % Set which timepoints ITC will be calculated over
foi    = [2 30];          % Set which frequencies ITC will be calculated over
chans  = 1:numel(ct_data.label); % Set channels to loop over (all is quite lengthy\0

% Swap time the time axis for brain time to avoid confusion for FieldTrip
% (cycles -> seconds)
for tr = 1:numel(bt_warpeddata.data.trial)
bt_warpeddata.data.time{tr} = ct_data.time{tr}(1:numel(bt_warpeddata.data.time{tr}));
end

% Use FieldTrip to calculate ITC
cfg           = [];
cfg.method    = 'wavelet';
cfg.toi       = toi;
cfg.output    = 'fourier';
cfg.foilim    = foi;

% Clock time
for c = chans % loop over channels
disp([num2str(round(c./numel(chans),2).*100),'% of Clock time ITC']);

cfg.channel = ct_data.label(c);
ct_dat      = ft_selectdata(cfg,ct_data);    
ct_freq     = ft_freqanalysis(cfg,ct_dat); 
ct_four     = ct_freq.fourierspctrm;        % copy the Fourier spectrum

ntrls     = size(ct_four,1);                % number of trials
ct_itc    = ct_four./abs(ct_four);          % divide by amplitude
ct_itc    = sum(ct_itc,1);                  % sum angles
ct_itc    = abs(ct_itc)/ntrls;              % take the absolute value and normalize
ct_itc    = squeeze(ct_itc);                % remove the first singleton dimension
ct_itc_full(c,:,:) = ct_itc;
end
ct_itc_mn = squeeze(mean(ct_itc_full,1));        % Average across channels

% Brain time
for c = chans % loop over channels
disp([num2str(round(c./numel(chans),2).*100),'% of Brain time ITC']);
    
cfg.channel = bt_warpeddata.data.label(c);
bt_dat      = ft_selectdata(cfg,bt_warpeddata.data);    
bt_freq     = ft_freqanalysis(cfg,bt_dat); 
bt_four     = bt_freq.fourierspctrm;        % copy the Fourier spectrum

ntrls     = size(bt_four,1);                % number of trials
bt_itc    = bt_four./abs(bt_four);          % divide by amplitude
bt_itc    = sum(bt_itc,1);                  % sum angles
bt_itc    = abs(bt_itc)/ntrls;              % take the absolute value and normalize
bt_itc    = squeeze(bt_itc);                % remove the first singleton dimension
bt_itc_full(c,:,:) = bt_itc;
end
bt_itc_mn = squeeze(mean(bt_itc_full,1));        % Average across channels

% Plot results
figure; hold on; subplot(2,1,1); 
imagesc(ct_freq.time, ct_freq.freq, ct_itc_mn);
axis xy
xlabel('time [s]'); ylabel('frequency [hz]');
title('Clock time ITC');
colorbar
hold on; subplot(2,1,2); 
imagesc(bt_freq.time, bt_freq.freq, bt_itc_mn);
axis xy
xlabel('time [s]'); ylabel('frequency [hz]');
title('Brain time ITC');
colorbar

% Depending on simulation paramers, the chosen warping signal,
% and selected channels we may see a stronger ITC in brain time than
% clock time.
% Results vary across parameters.

% To test the specificity of brain time warping, you may simulate a dataset
% and warp to frequencies absent in the data. Then, when repeating these
% scripts, the results should not show a significant difference between
% clock and brain time (unless a type I error occurs).