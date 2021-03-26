%%% Simulate a simple attention model, with three primary sources
%%% (dipoles). A conducter dipole in the right parietal cortex is
%%% conducting one follower dipole in each hemisphere's visual cortex.
%%% There are two conditions: attention to the left and right hemifield.
%%% The conductor dipole oscillates in phase with the follower dipole in
%%% the visual cortex contralateral to the attended hemifield, and in
%%% antiphase with the ipsilateral visual cortex dipole.

% Visualization of the simulation:
figure;imshow('bt_dipsim.png');

% establish parameters
settings.saveopt    = 1;          % Save the results?
settings.dur        = 3;          % Duration of simulated trials
settings.sr         = 200;        % Sampling rate of simulated data
settings.numtrials  = 60;         % Number of trial per condition
settings.mainfreq   = 10;         % Simulated frequency of interest
settings.freqdrift  = 0.05;       % Frequency walk per timepoint (random value between -n to +n)
settings.pinkns_amp = 0.5;        % Amplitude of 1/f pink noise
settings.savedir    = '\\its-rds.bham.ac.uk\rdsprojects\w\wimberm-ieeg-compute\Sander\TGM\BrainTimeToolbox-github\dipolesimulation';

% load relevant EEG files
load dip_elec
load dip_grid
load dip_vol

%% Generate signal parameters
t = linspace(0,settings.dur,settings.dur*settings.sr); % time vector
pinkns = pinknoise(numel(t),settings.pinkns_amp); % 1/f noise

% This simulation contains 11 dipoles - hardcoded so please don't change
num_active_dips = 11;

%% Generate random parameters for dipole generation to grab from
% Dipoles have:
% - Position [x y z]
% - Momentum [left-rightedness up-downness diagonality]
% - Frequency (x Hz)

% Generate random POSITION for 1/f dipoles
rand_pos_vec_x = randi([-55,60],num_active_dips,1);
rand_pos_vec_y = randi([-80,58],num_active_dips,1);
rand_pos_vec_z = randi([11,60],num_active_dips,1);
randpos = [rand_pos_vec_x rand_pos_vec_y rand_pos_vec_z];

% Generate random MOMENTUM for 1/f dipoles
rand_mom_vec = -1+(1+1)*rand(num_active_dips*3,1);
rand_mom_vec = reshape(rand_mom_vec,num_active_dips,3);
randmom = rand_mom_vec;

% Generate random FREQUENCIES for 1/f dipoles; but avoid frequency of interest
rand_freq_vec1 = randi([1 settings.mainfreq-2],1,floor(num_active_dips/2));
rand_freq_vec2 = randi([settings.mainfreq+2 20],1,ceil(num_active_dips/2));
randfreq = [rand_freq_vec1,rand_freq_vec2];
randfreq = randfreq(randperm(length(randfreq)));

% Generate random PHASE for 1/f dipoles
rand_phs=rand(settings.numtrials,num_active_dips)*2*pi;

% Generate random TRIAL START PHASE for conductor and follower dipoles
rand_phi=rand(settings.numtrials,num_active_dips)*2*pi;

% Generate FREQUENCY DRIFT for conductor and follower dipoles
% Employs a "random walk" type of approach
drift_sig = zeros(settings.numtrials,numel(t));

% Create drift
for trl = 1:settings.numtrials % Loop through trials
    for i = 1:numel(t)-1 % And timepoints
        r = -settings.freqdrift + (settings.freqdrift+settings.freqdrift)*rand(1,1); % Random drift
        drift_sig(trl,i+1) = drift_sig(trl,i)+r; % Add new drift to ongoing drift
    end
    drift_sig(trl,:) = smoothdata(drift_sig(trl,:),'Gaussian',settings.sr); % Smooth using 1-second window
end

% Plot drift for 6 example trials;
f1 = figure; set(gcf, 'WindowState', 'maximized');
subplot(5,5,[1 7]); hold on;
plot(t,drift_sig(1:6,:),'LineWidth',2);
title('Frequency drift over time (example trial 1 to 6)');
xlabel('Time (s)');
ylabel('Frequency change (Hz)');
set(gca,'FontName','Arial');
set(gca,'FontSize',16);


%% Main dipole signal creation
% Dipole 1 (CONDUCTOR dipole 1)
amp = 1;     % High amplitude
for tr=1:settings.numtrials
    dip_signal{tr}(1,:) = sin(2*pi*t.*(drift_sig(tr,:)+settings.mainfreq)+rand_phi(tr)+(pi/6)).*amp+pinkns; %+pi/6 to enable slight leading
end

% Dipole 2 (FOLLOWER dipole 1 [main])
amp = 1;     % High amplitude
offset = 0;  % Offset from zero?
for tr=1:settings.numtrials
    dip_signal{tr}(2,:) = (sin(2*pi*t.*(drift_sig(tr,:)+settings.mainfreq)+rand_phi(tr))+offset).*amp+pinkns;
end

% Dipole 3 (FOLLOWER dipole 2 [sub])
amp = 0.5;   % Medium amplitude
offset = 0;  % Offset from zero?
for tr=1:settings.numtrials
    dip_signal{tr}(3,:) = (sin(2*pi*t.*(drift_sig(tr,:)+settings.mainfreq)+rand_phi(tr)+pi)+offset).*amp+pinkns; %+pi yields antiphase relationship
end

% Plot primary dipoles
subplot(5,5,[16 22]),plot(t,[dip_signal{1}(1,:);dip_signal{1}(3,:);dip_signal{1}(2,:)],'LineWidth',2);
title('Signals of the three primary dipoles (example trial 1)');
xlabel('Time (s)');
ylabel('Amplitude (a.u.)');
legend('Conductor','Follower [main]','Follower [sub]');
set(gca,'FontName','Arial');
set(gca,'FontSize',16);

% 1/F pink noise dipoles
% Dipole 4 (1/F dipole 1)
for tr=1:settings.numtrials
    dip_signal{tr}(4,:) = cos(randfreq(1)*rand_phs(tr,1)*t)+pinkns; %random phase for this dipole
end

% Dipole 5 (1/F dipole 2)
for tr=1:settings.numtrials
    dip_signal{tr}(5,:) = sin(randfreq(2)*rand_phs(tr,2)*t)+pinkns;
end

% Dipole 6 (1/F dipole 3)
for tr=1:settings.numtrials
    dip_signal{tr}(6,:) = cos(randfreq(3)*rand_phs(tr,3)*t)+pinkns;
end

% Dipole 7 (1/F dipole 4)
for tr=1:settings.numtrials
    dip_signal{tr}(7,:) = sin(randfreq(4)*rand_phs(tr,4)*t)+pinkns;
end

% Dipole 8 (1/F dipole 5)
for tr=1:settings.numtrials
    dip_signal{tr}(8,:) = cos(randfreq(5)*rand_phs(tr,5)*t)+pinkns;
end

% Dipole 9 (1/F dipole 6)
for tr=1:settings.numtrials
    dip_signal{tr}(9,:) = sin(randfreq(6)*rand_phs(tr,6)*t)+pinkns;
end

% Dipole 10 (1/F dipole 7)
for tr=1:settings.numtrials
    dip_signal{tr}(10,:) = cos(randfreq(7)*rand_phs(tr,7)*t)+pinkns;
end

% Dipole 11 (1/F dipole 8)
for tr=1:settings.numtrials
    dip_signal{tr}(11,:) = sin(randfreq(8)*rand_phs(tr,8)*t)+pinkns;
end

%% Conducter and follower dipole parameters
% Look up brain area MNI locations: http://sprout022.sprout.yale.edu/mni2tal/mni2tal.html
% Right occipital cortex: [11 -92 12]
% Left occipital cortex: [-3 -93 12]
% Right parietal cortex: [46 -59 31] -- [50 -41 52] IPs
% Left parietal cortex: [-46 -60 33]

%% Let's start off with the "attended RIGHT" condition; where the left visual cortex is the main follower
% Dipole 1 (CONDUCTOR dipole 1)
dip1_pos = [50 -41 52]; %Right parietal
dip1_mom = [1 0 0]; % Momentum (left-rightedness, up-downness, diagonality)

% Dipole 2 (FOLLOWER dipole 1 [main])
dip2_pos = [-3 -93 12]; %Left occ
dip2_mom = [1 0 0]; % Momentum (left-rightedness, up-downness, diagonality)

% Dipole 2 (FOLLOWER dipole 2 [sub])
dip3_pos = [11 -92 12]; %Right occ
dip3_mom = [1 0 0]; % Momentum (left-rightedness, up-downness, diagonality)

% Put 1/F dipoles at random locations in the brain
% % Dipole 4 (1/F dipole 1)
dip4_pos = randpos(1,:);
dip4_mom = randmom(1,:);

% % Dipole 5 (1/F dipole 2)
dip5_pos = randpos(2,:);
dip5_mom = randmom(2,:);

% % Dipole 6 (1/F dipole 3)
dip6_pos = randpos(3,:);
dip6_mom = randmom(3,:);

% % Dipole 7 (1/F dipole 4)
dip7_pos = randpos(4,:);
dip7_mom = randmom(4,:);

% % Dipole 8 (1/F dipole 5)
dip8_pos = randpos(5,:);
dip8_mom = randmom(5,:);

% % Dipole 9 (1/F dipole 6)
dip9_pos = randpos(6,:);
dip9_mom = randmom(6,:);

% % Dipole 10 (1/F dipole 7)
dip10_pos = randpos(7,:);
dip10_mom = randmom(7,:);

% % Dipole 11 (1/F dipole 8)
dip11_pos = randpos(8,:);
dip11_mom = randmom(8,:);

%% Plot sensor locations
subplot(5,5,[3 25]); hold on;
p1 = ft_plot_mesh(dip_vol.bnd(1,3), 'facealpha', .1,'facecolor', 'skin');hold on;
p2 = ft_plot_mesh(dip_grid.pos(dip_grid.inside,:));
p3 = plot3(dip_elec.elecpos(:,1),dip_elec.elecpos(:,2),dip_elec.elecpos(:,3),'g*');hold on;

%% Create FieldTrip structure
for dip_i = 1:num_active_dips
    %Create position structure
    currdip_pos = ['dip',num2str(dip_i),'_pos'];
    dippos(dip_i,1:3) = eval(currdip_pos);
    dippos_struc{dip_i} = eval(currdip_pos);
    
    %Create momentum structure
    currdip_mom = ['dip',num2str(dip_i),'_mom'];
    dipmom(dip_i,1:3) = eval(currdip_mom);
    dipmom_struc{dip_i} = eval(currdip_mom);
    
    %Plot dipoles
    if dip_i <= 3
        p4 = plot3(dippos_struc{dip_i}(:,1),dippos_struc{dip_i}(:,2),dippos_struc{dip_i}(:,3),'ro','MarkerSize',20,'MarkerFaceColor','r');hold on;
    else
        p5 = plot3(dippos_struc{dip_i}(:,1),dippos_struc{dip_i}(:,2),dippos_struc{dip_i}(:,3),'bo','MarkerSize',10,'MarkerFaceColor','b');hold on;
    end
end

% Add plot information
legend([p3(1),p4(1),p5(1)],'Electrode','Primary Dipole','1/F Dipole');
set(gca,'FontName','Arial');
set(gca,'FontSize',16);

%% Simulate attended right data using FieldTrip
momvec = reshape(dipmom',1,numel(dipmom)); %Transform the momentum parameters into column vector

% Create FieldTrip structure
cfg=[];
cfg.dip.pos    = dippos;
cfg.dip.mom    = momvec';
cfg.relnoise   = 1;
cfg.vol        = dip_vol;
cfg.elec       = dip_elec;
cfg.dip.signal = dip_signal;
cfg.fsample    = settings.sr;
ct_right       = ft_dipolesimulation(cfg);

% Timelock the results
cfg=[];
cfg.keeptrials = 'yes';
cfg.covariance = 'yes';
cfg.covariancewindow = 'all';
ct_right=ft_timelockanalysis(cfg,ct_right);


%% Attended left condition; just switch around the two follower dipoles
% Dipole 2 (FOLLOWER dipole 1 [main])
dip2_pos = [11 -92 12]; %Right occ
dip2_mom = [1 0 0]; % Momentum (left-rightedness, up-downness, diagonality)

% Dipole 3 (FOLLOWER dipole 2 [sub])
dip3_pos = [-3 -93 12]; %Left occ
dip3_mom = [1 0 0]; % Momentum (left-rightedness, up-downness, diagonality)

for dip = 1:num_active_dips
    %Create position structure
    currdip_pos = ['dip',num2str(dip),'_pos'];
    dippos(dip,1:3) = eval(currdip_pos);
    dippos_struc{dip} = eval(currdip_pos);
    
    %Create momentum structure
    currdip_mom = ['dip',num2str(dip),'_mom'];
    dipmom(dip,1:3) = eval(currdip_mom);
    dipmom_struc{dip} = eval(currdip_mom);
end

%% Simulate attended left data using FieldTrip
% Create FieldTrip structure
cfg=[];
cfg.dip.pos    = dippos;
cfg.dip.mom    = momvec';
cfg.relnoise   = 1;
cfg.vol        = dip_vol;
cfg.elec       = dip_elec;
cfg.dip.signal = dip_signal;
cfg.fsample    = settings.sr;
ct_left       = ft_dipolesimulation(cfg);

% Timelock the results
cfg=[];
cfg.keeptrials = 'yes';
cfg.covariance = 'yes';
cfg.covariancewindow = 'all';
ct_left=ft_timelockanalysis(cfg,ct_left);

%% Now add condition-specific trialinfo that can be used as classification labels
ct_left.trialinfo = ones(size(ct_left.trial,1),1);
ct_right.trialinfo = 2*(ones(size(ct_right.trial,1),1));

% Save settings to trace back parameters
dipolesim_params = settings;

%% Save results
if settings.saveopt == 1
    save([settings.savedir,'\dipolesim'],'ct_left','ct_right','dipolesim_params','-v7.3');
    saveas(f1,[settings.savedir,'\dipolesim_positions.jpg'],'jpeg')
end
