function bt_wavplot(dat,ncycles,foi)


% Input:
% dat=NxM matrix containing the data to plot; rows=number of
% observations; columns = number of conditions
% col = color in rgb format (default: [0.8 0.2 0.5]);
% dotsize = size of dots to be plotted (default: 50)
% smooth = width of kernel used for smoothing (default: 100)

% Handle input for 1 frequency
if size(dat,2) == 1
    dat = dat';
end

% Transform tvec to phase
tvec = -ncycles/2:0.01:ncycles/2;
tvec = 2*pi.*tvec;

if ncycles~=2
    warning('X-axis will only be plotted in pi fractions for 2 cycles');
end

ncond=size(dat,1);
col = bt_colorscheme('symmetry');

x = linspace(foi(1)-0.6,foi(end)+0.6);
y = linspace(0,0);
plot(y,x,'linewidth',1,'color',[0.5 0.5 0.5],'linestyle', '--');
ylabel('Hz');
hold on

for k=1:numel(foi)
    tmp=dat(k,:);
    
    tmp_resc = rescale(tmp,foi(k)-0.4,foi(k)+0.4);
    plot(tvec,tmp_resc,'Color',col(k,:),'LineWidth',3);
    hold on
end

ylim([foi(1)-0.6 max(x)]);

ax=gca;
ax.YTick = foi;
ax.Box = 'off';
ax.XAxis.Exponent = 0;

if ncycles == 2
    xticks(-2*pi:pi:2*pi);
    xticklabels({'$-2\pi$','$-\pi$','0','$\pi$','$2\pi$'});
    set(gca,'TickLabelInterpreter', 'latex')
    xlabel('Phase (radians)');
else
    ylabel('Phase angle (degrees)');
end

% Adapt font
set(gca,'FontName',bt_plotparams('FontName'));
set(gca,'FontSize',bt_plotparams('FontSize'));
