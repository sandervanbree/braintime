function [comp] = bt_ica(config, data)

% Run ICA
cfg               = [];
cfg.method        = 'runica';
if isfield(config,'ncomps')
    cfg.runica.pca    = config.ncomps; % obtain N component, to reduce time
end
comp       = ft_componentanalysis(cfg ,data);


