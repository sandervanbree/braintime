function bt_figure(size)
% Color the figure background white and change its size (try two methods to handle old MATLAB versions).
% 'full' or no input = maximize figure
% 'half'             = intermediate size
% 'default'          = do not change size
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
    case 'half'
        set(gcf,'units','normalized','outerposition',[0.3 0.15 0.4 0.75]);
end

set(gcf,'color',bt_colorscheme('background'));