function bt_templatetopo(config)
% Create a template topography that will be used in bt_analyzecomps
% to bias component ranking toward the components' match to the desired
% topography. The template topography is automatically  saved in
% the /topography folder and loaded when using cfg.sortmethod = 'temptopo'
% in bt_analyzecomps.
%
% Use:
% bt_templatetopo(config)
%
% Input Arguments:
% config
%   - layout     % Layout file of the clock time data

% create a figure
f = figure;
cfg.layout  = config.layout;
layout = ft_prepare_layout(cfg);
ft_plot_lay(layout)

% add the required guidata
info       = guidata(gcf);
info.x     = layout.pos(:,1);
info.y     = layout.pos(:,2);
info.label = layout.label;
guidata(gcf, info)

    msgbox({'Draw a box around sensors of interest. Then, click in the box to continue.';' ';...
    'The component ranking will be biased by whether the component''s topography has an overlap with your sensors of interest';' ';...
    'The resulting template topography will be saved in the toolbox''s topography subfolder for retrieval when you enter cfg.sortmethod = ''templatetopo'' in bt_analyzecomps.'});

    % call to select_channel_ibm (slightly modified for the Brain Time Toolbox)
    set(gcf, 'WindowButtonDownFcn',   {@select_channel_ibm, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
    set(gcf, 'WindowButtonUpFcn',     {@select_channel_ibm, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
    set(gcf, 'WindowButtonMotionFcn', {@select_channel_ibm, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
    
    uiwait(f);
end