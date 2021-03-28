function bt_figure(maxopt)
% Color the figure background white and maximize (try two methods to handle old MATLABS)
% '1' or no input = maximize figure;
% '0' = do not maximize figure
% Color fetched from bt_colorscheme

switch nargin 
    case 0
        maxopt = 1; % Default
end

if maxopt == 1
    try set(gcf, 'WindowState', 'maximized');
    catch
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    end
end

set(gcf,'color',bt_colorscheme('background'));