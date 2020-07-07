function bt_templatetopo(config)

% create a figure
figure;
cfg.layout  = config.layout;
layout = ft_prepare_layout(cfg);
ft_plot_lay(layout)

% add the required guidata
info       = guidata(gcf);
info.x     = layout.pos(:,1);
info.y     = layout.pos(:,2);
info.label = layout.label;
guidata(gcf, info)

uiwait(msgbox({'Draw a box around sensors of interest. Then, click in the box to continue.';' ';...
    'The components will be ranked according to whether their topography has an overlap with your sensors of interest';' ';...
    'The resulting template topography will be saved in the toolbox''s template subfolder for later retrieval.'}));

    % call to select_channel_ibm (slightly modified for the Brain Time Toolbox)
    set(gcf, 'WindowButtonDownFcn',   {@select_channel_ibm, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
    set(gcf, 'WindowButtonUpFcn',     {@select_channel_ibm, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
    set(gcf, 'WindowButtonMotionFcn', {@select_channel_ibm, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
end