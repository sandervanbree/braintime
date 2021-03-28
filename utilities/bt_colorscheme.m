function [color] = bt_colorscheme(object)
% the Brain Time Toolbox fetches colors from this function, based on the
% specific object requested (input). You can manually set the color for
% different objects.
% In the first code block you can set line colors, in the second code block
% you can set color maps.
%
% Tip: Use 'uisetcolor' to browse colors and get their [r g b] values.

%% Lines and backgrounds (one color)
% Warping signal line indicator
warpingsignal = [0.38 0 0.45];

% Electrophysiology power spectra (not periodicity spectra)
ephys_ps = [0.45 0.45 0.45];

% Border lines of frequency range of interest
foi_borders = [0.3 0.3 0.3];

% Empirical periodicity spectra
per_ps_emp = [0.05 0.45 0.9];

% Permuted periodicity spectra
per_ps_perm = [0.4 0.4 0.4];

% Diagonal classification timecourse
diagonal = [0.04 0.23 0.03];

% P-value line in 1st and 2nd level statistics
pval = [1 0.68 0.12];

% Stars at significant frequencies in power spectra
% sigstar = [0.88 0.05 0.01];
sigstar = [1 0.68 0.12];

% Confidence interval fill color (is highly transparent)
confidenceinterval = [0.2 0.2 0.2]; 

% Significant cluster outline
sigcluster = [1 0.68 0.12];

% Non-significant cluster outline
nsigcluster = [0 0 0];

% Figure background
background = [1 1 1];

% Test color
% figure;plot(1:10,ones(1,10),'LineWidth',6,'Color',[0 1 0.2]);

%% Color maps (scale)
% Specify the code you would normally call in colormap

% Color map for topography plotting
topography = flipud(brewermap([],'RdBu'));
% topography = 'default';

% Color map for eletrophysiological time frequency plots
ephys_tfr = flipud(brewermap([],'RdBu'));

% Color map for Temporal Generalization Matrices
TGM = flipud(brewermap([],'RdBu'));

% Color map for Autocorrelation Maps
AC = flipud(brewermap([],'RdBu'));








%% Output requested map or color values
color = eval(object);
