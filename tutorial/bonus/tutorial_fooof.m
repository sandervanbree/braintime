%%% To determine the warping frequency of each participant, several methods
%%% of calculating individual power spectra may be used that account for
%%% aperiodic patterns on the data. A new method by Donoghue et al. (2020);
%%% Nat Neurosci. overcomes several limitations of classical approaches.
%%% This method has been implemented in the FOOOF toolbox for Python. In
%%% this tutorial, we will use FOOOF's MATLAB wrapper to extract peaks in
%%% power spectra from warping source data.
%%%
%%% To get this working, you need to have Python installed (versions later
%%% than 3.7 may not work) https://www.python.org/downloads/release/python-370/
%%% You also need to pip install fooof
%%% For all instructions, please go to the README file in the
%%% ...\external\fooof_mat\ folder.

% Step 1: Check whether Python is installed and accessible by MALAB
pyversion

% Load your warping source structure
load tutorial1_output warpsources 

% Restrict tutorial to 10 trials to speed things up
cfg         = [];
cfg.trials  = 1:10;
warpsources = ft_selectdata(cfg,warpsources);

% Get basic parameters
sr       = warpsources.fsample;                % Sampling rate
len      = numel(warpsources.trial{1}(1,:));   % Trial length
nwarpsrc = size(warpsources.trial{1},1);       % Number of warping sources
ntrls    = numel(warpsources.trial);           % Number of trials
cnt      = 1;                                  % Initiate counter 

% FOOOF settings
settings = struct();  % Use defaults
f_range  = [2, 30];   % Specify frequency range of interest

% FOOOF requires the data to be formatted time x channel

for trl = 1:ntrls                                           % Loop over all trial every trial
    
    disp(['FOOOFing trial ',num2str(trl),...               % Notify user of progress
        ' out of ',num2str(numel(warpsources.trial))]); 
    
    for ws = 1:nwarpsrc                                     % Loop over all warping sources
        currtrl = warpsources.trial{trl}(ws,:)';
        [psd, freqs] = pwelch(currtrl, len, [], [], sr);    % Get power spectrum                                                 
        freqs = freqs';                                     % Transpose, to make inputs row vectors
        psd = psd';
        
        % Run FOOOF, also returning the model
        fooof_results = fooof(freqs, psd, f_range, settings, true);
        
        ps(cnt,:)     = fooof_results.power_spectrum;       % Extract power spectrum
        model(cnt,:)  = fooof_results.fooofed_spectrum;     % Extract FOOOF model
        ap_fit(cnt,:) = fooof_results.ap_fit;               % Extract aperiodic fit
        freq(cnt,:)   = fooof_results.freqs;                % Extract frequency vector
        
        cnt = cnt+1;
    end
end

% Average all results
ps_avg       = mean(ps,1);
model_avg    = mean(model,1);
ap_fit_avg   = mean(ap_fit,1);
freq_avg     = mean(freq,1);  % All frequency vectors should be the same

% Plot the resulting model (code modified from fooof_plot.m)
figure;hold on;title('FOOOF results');
p = plot(freq_avg,ps_avg,'LineWidth',2,'Color','k');
m = plot(freq_avg,model_avg,'LineWidth',2,'Color','r');
a = plot(freq_avg,ap_fit_avg,'LineWidth',2,'Color','b','LineStyle','--');
m.Color(4) = 0.5;
grid on
xlabel('Frequency (Hz)'); ylabel('Power');
legend('Original Spectrum', 'Full Model Fit', 'Aperiodic Fit');

% Find peak in model spectrum
[~,pk_inds] = findpeaks(model_avg);

% Filter peaks to frequencies of interest
foi = 8:12;                                                        % Range of interest
foi_inds = find(freq_avg>=foi(1) & freq_avg<=foi(end));            % Find indices of range
pk_inds = pk_inds(pk_inds>=foi_inds(1) & pk_inds<=foi_inds(end));  % Filter peaks to range
pk_freqs = freq_avg(pk_inds);                                      % Get frequencies

% Plot potential warping frequencies
for i = numel(pk_freqs)
xline(pk_freqs(i),'LineWidth',2,'Color','g');
end
legend('Original Spectrum', 'Full Model Fit', 'Aperiodic Fit','Potential warping frequency');