function bt_rainplot(dat,col,dotsize,smooth,xlimit)
% Written by Simon Hanslmayr. Modified for the Brain Time Toolbox by Sander
% van Bree

% Input:
% dat=NxM matrix containing the data to plot; rows=number of
% observations; columns = number of conditions
% col = color in rgb format (default: [0.8 0.2 0.5]);
% dotsize = size of dots to be plotted (default: 50)
% smooth = width of kernel used for smoothing (default: 100)

if isempty(col)
    col = bt_colorscheme('asymm_wavshap');
end

if isempty(dotsize)
    dotsize = 50;
end

if isempty(smooth)
    smooth=100;
end



x = linspace(0-0.6,0+0.6);
y = linspace(0,0);
plot(y,x,'linewidth',1,'color',[0.5 0.5 0.5],'linestyle', '--');
ylabel('Hz');
xlabel('Asymmetry index');
hold on

    tmp=dat;
    jitt=(rand(1,length(tmp))+0.1).*0.1;
    scatter(tmp',jitt,dotsize,col,'filled','MarkerFaceAlpha',.5);
    hold on


ylim([0-0.6 max(x)]);


% Plot Boxes and means ...
xb=[0-0.1350 0-0.02];
xbr=[0-0.1350 0-0.02];


    Mns=mean(dat);
    STE=std(dat)./sqrt(length(dat));%
    STD=std(dat);%
    
    % draw box for STD
    tmpxb=[xb(1,1) xb(1,1) xb(1,2) xb(1,2) xb(1,1)];
    tmpyb=[Mns-STD Mns+STD Mns+STD Mns-STD Mns-STD];
    plot(tmpyb,tmpxb,'linewidth',1,'color',col,'linestyle', '-','LineWidth',1.5);
    % draw filled box for STE
    tmpxb=[xb(1,1) xb(1,1) xb(1,2) xb(1,2)];
    tmpyb=[Mns-STE Mns+STE Mns+STE Mns-STE];
    fill(tmpyb,tmpxb,col,'FaceAlpha',.8,'EdgeAlpha',0);
    hold on
    
    % Plot solid Lines for means
    x = linspace(xb(1,1),xb(1,2))-1;
    y = linspace(mean(Mns),mean(Mns));
    
    plot(y,x,'linewidth',1,'color',col,'linestyle', '-','LineWidth',2);
    hold on
    
    % Plot Error bars for 90th percentile
    pct = prctile(dat,[5 95]);
    x = linspace(xbr(1,1),xbr(1,2));
    y = linspace(pct(1),pct(1));
    % draw lower bound
    plot(y,x,'linewidth',1,'color',col,'linestyle', '-','LineWidth',1.5);
    hold on
    
    x = linspace(xbr(1,1),xbr(1,2));
    y = linspace(pct(2),pct(2));
    % draw upper bound
    plot(y,x,'linewidth',1,'color',col,'linestyle', '-','LineWidth',1.5);
    hold on
    
    % draw connecting lines to upper bound
    x = [mean(xbr,2),mean(xbr,2)];
    y = [Mns+STD,pct(2)];
    % draw upper bound
    plot(y,x,'linewidth',1,'color',col,'linestyle', '-','LineWidth',1.5);
    hold on
    
    % draw connecting lines to lower bound
    x = [mean(xbr,2),mean(xbr,2)];
    y = [Mns-STD,pct(1)];
    % draw upper bound
    plot(y,x,'linewidth',1,'color',col,'linestyle', '-','LineWidth',1.5);
    hold on
    


% Let it rain ...
nbins=400;
scal=0.3;
xpos=xb(1,1);
%xb=[0.2 0.2+scal;1.2 1.2+scal;2.2 2.2+scal];

    histdat=hist(dat,nbins);
    dens=conv(histdat,gausswin(smooth),'same')./(smooth/1000);% compute smoothed firing rate using gaussian kernel
    %dens=dens-min(dens);
    dens=(dens./(max(dens))).*-scal;
    tmpxb=[0 dens 0 0]+xpos*1.1;
    tmpyb=[min(dat) linspace(min(dat),max(dat),nbins) max(dat) min(dat)];
    %figure;plot(tmpxb(1:end),tmpyb(1:end));
    fill(tmpyb,tmpxb,col,'FaceAlpha',.8,'EdgeAlpha',0);
    hold on

ax=gca;
ax.YTick = 1;
ax.Box = 'off';
ax.XLim = xlimit;
ax.YAxis.Color = [1 1 1];
ax.XAxis.Exponent = 0;
set(gca,'TickLabelInterpreter', 'latex')
% Adapt font
set(gca,'FontName',bt_plotparams('FontName'));
set(gca,'FontSize',bt_plotparams('FontSize'));
