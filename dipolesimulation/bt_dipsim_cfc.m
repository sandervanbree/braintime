%%% Simulate a Cross-Frequency Coupling (CFC) model with two coupled
%%% oscillators in the left hemisphere, and two in the right. For each of
%%% the two conditions, the left two nested oscillators are in-phase and
%%% the right oscillators are in-phase; but the relation between the left
%%% and right set flips by 180 degrees as a function of condition. This
%%% enables us to test to what extent brain time warping can tap into
%%% nested oscillations (i.e., does warping to one oscillator bring out the
%%% other?).
%%%

% Visualization of the simulation:
figure;imshow('bt_dipsim_cfc.png');

% establish parameters
settings.saveopt    = 1;         % Save the results?
settings.dur        = 3;         % Duration of simulated trials
settings.sr         = 200;      % Sampling rate of simulated data
settings.numtrials  = 60;        % Number of trial per condition
settings.lowfreq    = 10;        % Frequency of the low freq oscillation
settings.highfreq   = 17;        % Frequency of the high freq nested oscillation
settings.freqdrift  = 0.05;      % 0.05 Frequency walk per timepoint (random value between -n to +n)
settings.pinkns_amp = 0.5;       % 0.5 Amplitude of 1/f pink noise
settings.randphase  = 1;         % Randomise starting phase of primary dipoles per trial?
settings.cfcphase   = 1;         % Cross-frequency PHASE coupling in additiont to Cross-frequency AMPLITUDE coupling?
                                 % for details: check Jensen & Colgin (2007); TICS
settings.relnoise   = 2.5;         % How much noise to add in ft_dipolesimulation
settings.savedir    = '\\analyse4\project0318\Sander\software\braintime-master\dipolesimulation';

% load relevant EEG files
load dip_elec
load dip_grid
load dip_vol

%% Generate signal parameters
t = linspace(0,settings.dur,settings.dur*settings.sr); % time vector
pinkns = pinknoise(numel(t),1); % 1/f noise
pinkns = pinkns.*settings.pinkns_amp; % Change 1/f amplitude

% This simulation contains 10 dipoles - hardcoded so please don't change
num_active_dips = 10;

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
rand_freq_vec1 = randi([1 settings.lowfreq-2],1,floor(num_active_dips/2));
rand_freq_vec2 = randi([settings.lowfreq+2 20],1,ceil(num_active_dips/2));
randfreq = [rand_freq_vec1,rand_freq_vec2];
randfreq = randfreq(randperm(length(randfreq)));

% Generate random PHASE for 1/f dipoles
rand_phs=rand(settings.numtrials,num_active_dips)*2*pi;

if settings.randphase == 1
% Generate random TRIAL START PHASE for conductor and follower dipoles
rand_phi=rand(settings.numtrials,num_active_dips)*2*pi;
else
% Generate stable TRIAL START PHASE for conductor and follower dipoles
rand_phi=ones(settings.numtrials,num_active_dips)*2*pi;
end

% Generate FREQUENCY DRIFT for conductor and follower dipoles
% Employs a "random walk" type of approach
drift_low = zeros(settings.numtrials,numel(t));
drift_high = zeros(settings.numtrials,numel(t));

% Create drift LOW FREQ oscillator
for trl = 1:settings.numtrials % Loop through trials
    for i = 1:numel(t)-1 % And timepoints
        r = -settings.freqdrift + (settings.freqdrift+settings.freqdrift)*rand(1,1); % Random drift
        drift_low(trl,i+1) = drift_low(trl,i)+r; % Add new drift to ongoing drift
    end
    drift_low(trl,:) = smoothdata(drift_low(trl,:),'Gaussian',settings.sr); % Smooth using 1-second window
end

% Create drift HIGH FREQ oscillator
for trl = 1:settings.numtrials % Loop through trials
    for i = 1:numel(t)-1 % And timepoints
        r = -settings.freqdrift + (settings.freqdrift+settings.freqdrift)*rand(1,1); % Random drift
        drift_high(trl,i+1) = drift_high(trl,i)+r; % Add new drift to ongoing drift
    end
    drift_high(trl,:) = smoothdata(drift_high(trl,:),'Gaussian',settings.sr); % Smooth using 1-second window
end

%% Main dipole signal creation
% Dipole 1 (RIGHT lowfreq dipole)
amp = 1;     % High amplitude
offset = 0;  % Offset from zero?
for tr=1:settings.numtrials
    dip_signal{tr}(1,:) = (sin(2*pi*t.*(drift_low(tr,:)+settings.lowfreq)+rand_phi(tr))+offset).*amp+pinkns;
end

% Dipole 2 (LEFT lowfreq dipole)
amp = 0.5;   % Medium amplitude
offset = 0;  % Offset from zero?
for tr=1:settings.numtrials
    dip_signal{tr}(2,:) = (sin(2*pi*t.*(drift_low(tr,:)+settings.lowfreq)+rand_phi(tr)+pi)+offset).*amp+pinkns; %+pi yields antiphase relationship
end

% Dipole 3 (RIGHT highfreq dipole)
amp = 1;     % High amplitude
offset = 0;  % Offset from zero?
for tr=1:settings.numtrials
    % If you include CFC phase coupling, this applies the same frequency
    % drift to BOTH nested oscillators.
    if settings.cfcphase == 1
    dip_signal{tr}(3,:) = (sin(2*pi*t.*(drift_low(tr,:)+settings.highfreq)+rand_phi(tr))+offset).*amp+pinkns;
    else
    dip_signal{tr}(3,:) = (sin(2*pi*t.*(drift_high(tr,:)+settings.highfreq)+rand_phi(tr))+offset).*amp+pinkns;    
    end
    % Now apply phase-amplitude coupling (via element-wise multiplication)
    dip_signal{tr}(3,:) = dip_signal{tr}(3,:).*(rescale((dip_signal{tr}(1,:)-squeeze(mean(dip_signal{tr}(1,:))))));
end

% Dipole 4 (LEFT highfreq dipole)
amp = 0.5;   % Medium amplitude
offset = 0;  % Offset from zero?
for tr=1:settings.numtrials
    if settings.cfcphase == 1
    dip_signal{tr}(4,:) = (sin(2*pi*t.*(drift_low(tr,:)+settings.highfreq)+rand_phi(tr)+pi)+offset).*amp+pinkns;
    else
    dip_signal{tr}(4,:) = (sin(2*pi*t.*(drift_high(tr,:)+settings.highfreq)+rand_phi(tr)+pi)+offset).*amp+pinkns;    
    end
    dip_signal{tr}(4,:) = dip_signal{tr}(4,:).*(rescale(dip_signal{tr}(2,:)-squeeze(mean(dip_signal{tr}(2,:)))));
end

% Plot primary dipoles
figure;subplot(2,1,1);plot(t,[dip_signal{1}(1,:);dip_signal{1}(3,:);],'LineWidth',2);
title('Right nested oscillators');
legend('Low freq','High freq');
xlabel('Time (s)');
ylabel('Amplitude (a.u.)');
set(gca,'FontName','Arial');
set(gca,'FontSize',16);
hold on;subplot(2,1,2);plot(t,[dip_signal{1}(2,:);dip_signal{1}(4,:)],'LineWidth',2);
title('Left nested oscillators');
legend('Low freq','High freq');
xlabel('Time (s)');
ylabel('Amplitude (a.u.)');
set(gca,'FontName','Arial');
set(gca,'FontSize',16);
 
% 1/F pink noise dipoles
% Dipole 5 (1/F dipole 1)
for tr=1:settings.numtrials
    dip_signal{tr}(5,:) = cos(randfreq(3)*rand_phs(tr,3)*t)+pinkns;
end

% Dipole 6 (1/F dipole 2)
for tr=1:settings.numtrials
    dip_signal{tr}(6,:) = sin(randfreq(4)*rand_phs(tr,4)*t)+pinkns;
end

% Dipole 7 (1/F dipole 3)
for tr=1:settings.numtrials
    dip_signal{tr}(7,:) = cos(randfreq(5)*rand_phs(tr,5)*t)+pinkns;
end

% Dipole 8 (1/F dipole 4)
for tr=1:settings.numtrials
    dip_signal{tr}(8,:) = sin(randfreq(6)*rand_phs(tr,6)*t)+pinkns;
end

% Dipole 9 (1/F dipole 5)
for tr=1:settings.numtrials
    dip_signal{tr}(9,:) = cos(randfreq(7)*rand_phs(tr,7)*t)+pinkns;
end

% Dipole 10 (1/F dipole 6)
for tr=1:settings.numtrials
    dip_signal{tr}(10,:) = sin(randfreq(8)*rand_phs(tr,8)*t)+pinkns;
end


%% Let's start off with one condition: RIGHT
% let's call it "right" because the right hemisphere oscillators are at
% 0pi, whereas the left hemisphere is at +pi

% Dipole 1 (RIGHT lowfreq dipole)
dip1_pos = [50 -41 52]; % RIGHT hemisphere
dip1_mom = [1 0 0]; % Momentum (left-rightedness, up-downness, diagonality)

% Dipole 2 (LEFT lowfreq dipole)
dip2_pos = [-50 -41 52]; % LEFT hemisphere
dip2_mom = [1 0 0]; % Momentum (left-rightedness, up-downness, diagonality)

% % Dipole 3 (RIGHT highfreq dipole)
dip3_pos = dip1_pos;
dip3_pos(2) = dip3_pos(2)+45; % keep it close (push it 45 units more rostral)
dip3_mom = dip1_mom;

% % Dipole 4 (LEFT highfreq dipole)
dip4_pos = dip2_pos;
dip4_pos(2) = dip4_pos(2)+45;
dip4_mom = dip2_mom;

% Put 1/F dipoles at random locations in the brain
% % Dipole 5 (1/F dipole 1)
dip5_pos = randpos(3,:);
dip5_mom = randmom(3,:);

% % Dipole 6 (1/F dipole 2)
dip6_pos = randpos(4,:);
dip6_mom = randmom(4,:);

% % Dipole 7 (1/F dipole 3)
dip7_pos = randpos(5,:);
dip7_mom = randmom(5,:);

% % Dipole 8 (1/F dipole 4)
dip8_pos = randpos(6,:);
dip8_mom = randmom(6,:);

% % Dipole 9 (1/F dipole 5)
dip9_pos = randpos(7,:);
dip9_mom = randmom(7,:);

% % Dipole 10 (1/F dipole 6)
dip10_pos = randpos(8,:);
dip10_mom = randmom(8,:);

%% Plot sensor locations
subplot(5,5,[3 25]); hold on;
p1 = ft_plot_mesh(dip_vol.bnd(1,3), 'facealpha', .1,'facecolor', 'skin');hold on;
p2 = ft_plot_mesh(dip_grid.pos(dip_grid.inside,:));
% p3 = plot3(dip_elec.elecpos(:,1),dip_elec.elecpos(:,2),dip_elec.elecpos(:,3),dip_elec.elecpos(:,4),dip_elec.elecpos(:,5),'g*');hold on;

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
    if dip_i < 3
        p4 = plot3(dippos_struc{dip_i}(:,1),dippos_struc{dip_i}(:,2),dippos_struc{dip_i}(:,3),'ro','MarkerSize',20,'MarkerFaceColor','r');hold on;
    elseif dip_i > 2 && dip_i < 5
        p5 = plot3(dippos_struc{dip_i}(:,1),dippos_struc{dip_i}(:,2),dippos_struc{dip_i}(:,3),'ro','MarkerSize',20,'MarkerFaceColor','magenta');hold on;        
    else
        p6 = plot3(dippos_struc{dip_i}(:,1),dippos_struc{dip_i}(:,2),dippos_struc{dip_i}(:,3),'bo','MarkerSize',10,'MarkerFaceColor','b');hold on;
    end
end

% Add plot information
legend([p4(1),p5(1),p6(1)],'low freq Dipole','high freq dipole','1/F Dipole');
set(gca,'FontName','Arial');
set(gca,'FontSize',16);

%% Simulate right class data using FieldTrip
momvec = reshape(dipmom',1,numel(dipmom)); %Transform the momentum parameters into column vector

% Create FieldTrip structure
cfg=[];
cfg.dip.pos    = dippos;
cfg.dip.mom    = momvec';
cfg.relnoise   = settings.relnoise;
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


%% Left condition; just switch around the dipoles (now the right hemisphere will have +pi)
% Dipole 1 (RIGHT lowfreq dipole)
dip1_pos = [-50 -41 52]; % RIGHT hemisphere
dip1_mom = [1 0 0]; % Momentum (left-rightedness, up-downness, diagonality)

% Dipole 2 (LEFT lowfreq dipole)
dip2_pos = [50 -41 52]; % LEFT hemisphere
dip2_mom = [1 0 0]; % Momentum (left-rightedness, up-downness, diagonality)

% % Dipole 3 (RIGHT highfreq dipole)
dip3_pos = dip1_pos;
dip3_pos(2) = dip3_pos(2)+45; % keep it close (push it 45 units more rostral)
dip3_mom = dip1_mom;

% % Dipole 4 (LEFT highfreq dipole)
dip4_pos = dip2_pos;
dip4_pos(2) = dip4_pos(2)+45;
dip4_mom = dip2_mom;

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
cfg.relnoise   = settings.relnoise;
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
dipolesim_params          = settings;
dipolesim_params.dipfreqs = randfreq(3:8); % add simulated 1/f frequencies

% Save data
if settings.saveopt == 1
    save([settings.savedir,'\dipolesim_cfc'],'ct_left','ct_right','dipolesim_params','-v7.3');
end

