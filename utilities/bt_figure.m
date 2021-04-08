function bt_figure(size)
% Color the figure background white and change its size (try two methods to handle old MATLAB versions).
% 'full' or no input = maximize figure
% 'default'          = do not change size
% 'halflong'         = intermediate size (vertical orientation)
% 'halfwide'         = intermediate size (horizontal orientation)
% 'cluster'          = suitable size for cluster output
%
% Color fetched from bt_colorscheme

switch nargin
    case 0
        size = 'full'; % Default
end

switch size
    case 'full'
        try set(gcf, 'WindowState', 'maximized');
        catch
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        end
    case 'halflong'
        set(gcf,'units','normalized','outerposition',[0.3 0.15 0.4 0.75]);
    case 'halfwide'
        set(gcf,'units','normalized','outerposition',[0.2 0.15 0.55 0.75]);
    case 'cluster'
        set(gcf,'units','normalized','outerposition',[0.2 0.15 0.55 0.75])
    case 'braintime_per'
        set(gcf,'units','normalized','outerposition',[0.0948 0.1630 0.8047 0.6639])
    case 'clocktime_per'
        set(gcf,'units','normalized','outerposition',[0.2 0.15 0.55 0.75]);
end

set(gcf,'color',bt_colorscheme('background'));