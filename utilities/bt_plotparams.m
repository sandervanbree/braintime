function [param] = bt_plotparams(object)
% the Brain Time Toolbox fetches plotting parameters from this function,
% based on the specific object requested (input). You can manually set the
% parameters for different objects.

%% Font
% Font (except for bt_selectsource)
FontName = 'Arial';

% Font size (except for bt_selectsource)
FontSize = 15.5;

%% Axis limits 
% By default, temporal generalization matrices (TGM), autocorrelation maps
% (AC), and diagonal classification timecourses have a data-driven limit.
% You can override this by entering a minimum and maximum (e.g. 
% [0.35,0.65]). This may come in handy when visually comparing different
% participants.
% Leave an object empty ([]) to enable data-driven method.

% Color lmits of TGM
TGM_clim = []; % [0.35 0.65]

% Color limits of AC
AC_clim = [];

% Y-axis limits of diagonal
diag_ylim = [];

%% Fetch value from requested object
param = eval(object);