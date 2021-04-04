function bt_rainplot_legacy(dat,col,dotsize,smooth,foi)
% This is a legacy function that is no longer part of the Brain Time
% Toolbox, but may come in useful.
% Please see tutorial_checksymmetry and bt_checksymmetry for details.
%
% Written by Simon Hanslmayr. Modified for the Brain Time Toolbox by Sander
% van Bree
%
% Input:
% dat=NxM matrix containing the data to plot; rows=number of
% observations; columns = number of conditions
% col = color in rgb format (default: [0.8 0.2 0.5]);
% dotsize = size of dots to be plotted (default: 50)
% smooth = width of kernel used for smoothing (default: 100)

ncond=size(dat,1);
if isempty(col)
    col = bt_colorscheme('checksymmetry');
end

if isempty(dotsize)
    dotsize = 50;
end

if isempty(smooth)
    smooth=100;
end

x = linspace(foi(1)-0.6,foi(end)+0.6);
y = linspace(0,0);
plot(y,x,'linewidth',1,'color',[0.5 0.5 0.5],'linestyle', '--');
ylabel('Hz');
xlabel('Asymmetry index');
hold on

for k=1:numel(foi)
    tmp=dat(k,:);
    jitt=(rand(1,length(tmp))+0.2).*0.2;
    scatter(tmp,foi(k)+jitt,dotsize,col(k,:),'filled','MarkerFaceAlpha',.6);
    hold on
end

ylim([foi(1)-0.6 max(x)]);


% Plot Boxes and means ...
xb=[foi(1)-0.1350 foi(1)-0.02];
xbr=[foi(1)-0.1350 foi(1)-0.02];

for plt=1:k
    Mns=mean(dat(plt,:));
    STE=std(dat(plt,:))./sqrt(length(dat));%
    STD=std(dat(plt,:));%
    
    % draw box for STD
    tmpxb=[xb(1,1) xb(1,1) xb(1,2) xb(1,2) xb(1,1)]+plt-1;
    tmpyb=[Mns-STD Mns+STD Mns+STD Mns-STD Mns-STD];
    plot(tmpyb,tmpxb,'linewidth',1,'color',col(plt,:),'linestyle', '-','LineWidth',1.5);
    % draw filled box for STE
    tmpxb=[xb(1,1) xb(1,1) xb(1,2) xb(1,2)]+plt-1;
    tmpyb=[Mns-STE Mns+STE Mns+STE Mns-STE];
    fill(tmpyb,tmpxb,col(plt,:),'FaceAlpha',.6,'EdgeAlpha',0);
    hold on
    
    % Plot solid Lines for means
    x = linspace(xb(1,1)+plt,xb(1,2)+plt)-1;
    y = linspace(mean(Mns),mean(Mns));
    
    plot(y,x,'linewidth',1,'color',col(plt,:),'linestyle', '-','LineWidth',2);
    hold on
    
    % Plot Error bars for 90th percentile
    pct = prctile(dat(plt,:),[5 95]);
    x = linspace(xbr(1,1)+plt,xbr(1,2)+plt)-1;
    y = linspace(pct(1),pct(1));
    % draw lower bound
    plot(y,x,'linewidth',1,'color',col(plt,:),'linestyle', '-','LineWidth',1.5);
    hold on
    
    x = linspace(xbr(1,1)+plt,xbr(1,2)+plt)-1;
    y = linspace(pct(2),pct(2));
    % draw upper bound
    plot(y,x,'linewidth',1,'color',col(plt,:),'linestyle', '-','LineWidth',1.5);
    hold on
    
    % draw connecting lines to upper bound
    x = [mean(xbr,2),mean(xbr,2)]+plt-1;
    y = [Mns+STD,pct(2)];
    % draw upper bound
    plot(y,x,'linewidth',1,'color',col(plt,:),'linestyle', '-','LineWidth',1.5);
    hold on
    
    % draw connecting lines to lower bound
    x = [mean(xbr,2),mean(xbr,2)]+plt-1;
    y = [Mns-STD,pct(1)];
    % draw upper bound
    plot(y,x,'linewidth',1,'color',col(plt,:),'linestyle', '-','LineWidth',1.5);
    hold on
    
end

% Let it rain ...
nbins=400;
scal=0.3;
xpos=0.5;
%xb=[0.2 0.2+scal;1.2 1.2+scal;2.2 2.2+scal];

for plt=1:ncond
    histdat=hist(dat(plt,:),nbins);
    dens=conv(histdat,gausswin(smooth),'same')./(smooth/1000);% compute smoothed firing rate using gaussian kernel
    %dens=dens-min(dens);
    dens=(dens./(max(dens))).*-scal;
    tmpxb=[0 dens 0 0]-min(dens)+foi(plt)-1+xpos;
    tmpyb=[min(dat(plt,:)) linspace(min(dat(plt,:)),max(dat(plt,:)),nbins) max(dat(plt,:)) min(dat(plt,:))];
    %figure;plot(tmpxb(1:end),tmpyb(1:end));
    fill(tmpyb,tmpxb,col(plt,:),'FaceAlpha',.6,'EdgeAlpha',0);
    hold on
end

ax=gca;
ax.YTick = foi;
ax.Box = 'off';
ax.XAxis.Exponent = 0;
set(gca,'TickLabelInterpreter', 'latex')
% Adapt font
set(gca,'FontName',bt_plotparams('FontName'));
set(gca,'FontSize',bt_plotparams('FontSize'));
