function [bt_source] = bt_selectsource(config, fft_sources, warpsources)
% Display the warping sources analyzed and ranked by bt_analyzesources.
% This allows the user to visually inspect warping sources. Each warping
% source contains a potential warping signal which can be used to brain
% time warp the data. Please select the optimal warping source based on the
% time-frequency characteristics of the warping signal, the topography of
% the warping source (if relevant), and several other criteria described
% in the function's graphical interface.
%
% Use:
% [bt_source] = bt_selectsource(cfg,fft_sources,warpsources)
%
% Input Arguments:
% config
%   - layout         % (Optional) A FieldTrip layout file to plot warping
%                    % source data at the sensor level. See
%                    % ft_prepare_layout for details.
%                    %
% fft_sources        % Data structure with time frequency information of
%                    % all warping sources as obtained by
%                    % bt_analyzesources.
%                    %
% warpsources        % FieldTrip data structure that contains warping
%                    % sources (ICA components, virtual channels,
%                    % or LFP time). Each warping source contains a
%                    % warping signal, one of which will be selected here.
%                    %
% Output:            %
% bt_source          % Data structure with: selected warping signal, its
%                    % time frequency information, and config details

%% Get basic info
srcrank = fft_sources{1};                              % Time freq data of the ranked sources
mintime = fft_sources{2}(1);                           % Start time of interest
maxtime = fft_sources{2}(2);                           % End time of interest
minfoi = fft_sources{3}(1);                            % Lowest possible warping frequency
maxfoi = fft_sources{3}(2);                            % Highest possible warping frequency
fspecinfo = fft_sources{4};                            % FFT time and frequency vector
minfft = fspecinfo.freq(1);                            % Lowest analyzed freq in FFT
maxfft = fspecinfo.freq(end);                          % Highest analyzed freq in FFT
powtf = fft_sources{5};                                % Power spectrum of warping sources
pspec = fft_sources{6};                                % Power spectrum averaged across trials
phs = fft_sources{7};                                  % Phase of all warping sources for all trials
cutmethod = fft_sources{8};                            % Applied cutting method
rankmethod = fft_sources{9};                           % Applied ranking method

mintime_ind = nearest(fspecinfo.time,mintime);             % Index of start time of interest (differs for cutartefact)
maxtime_ind = nearest(fspecinfo.time,maxtime);             % Index of end time of interest


%% Plot ranked warping sources
% Find out ranking type
if strcmp(rankmethod,'maxpow')
    ranktype = 'maximum power.';
elseif strcmp(rankmethod,'templatetopo')
    ranktype = 'maximum power and correlation to template topography.';
end

% Prepare instruction message
msg = ['The warping sources are ranked by ',ranktype,newline,...
    'Please select the warping source containing your warping signal, based on its:',newline,...
    '(1) time frequency characteristics',newline,...
    '(2) topography of interest (if applicable)',newline,...
    '(3) variance explained (for ICA components, lower source number means higher r^2)',newline,newline,...
    'INSTRUCTIONS:',newline,...
    'Press \leftarrow or \rightarrow arrow to see the previous or next warping source. Click on a',newline,...
    'warping source to highlight it, then keep browsing or press ''Q'', spacebar, or enter to quit.'];

% Prepare figure
src_oi = -2; %source of interest
maxp=max(max(max(powtf(:,:,:)))); %establish maximum power. This will be used to limit axes and normalize power to
% Sanity check
if isnan(maxp)
    error('The time frequency analysis of the warping source analysis has yielded only NaNs. Did you select the right time window?, and does your data contain additional time before and after your time window?');
end
finish = 0;
src_ind=1;
f1 = figure;hold on;

% Adapt figure
bt_figure;

% Create loop to browse through warping sources
while finish==0
    if src_ind < 1 % make sure source number cannot go out of bounds
        src_ind = 1;
    elseif src_ind>size(srcrank,1)
        src_ind = size(srcrank,1);
    end
    
    % Establish current warping source and its maximum frequency
    currsrc = srcrank(src_ind);
    maxfreq = srcrank(src_ind,2);
    
    subplot(8,2,[1 3]);
    inst = text(0,0.5,msg,'FontName','Arial','FontSize',11.5); axis off
    
    sp1 = subplot(8,2,[7 9 11 13 15]);
    % Only plot topography if layout is specified
    if isfield(config,'layout')
        % warping source topography
        cfg           = [];
        cfg.figure    = sp1;
        cfg.component = currsrc; % specify the source(s) that should be plotted
        cfg.layout    = config.layout; % specify the layout file that should be used for plotting
        cfg.comment   = 'no';
        
        % "warping source" could be either ICA comps or another data structure
        try
            ft_topoplotIC(cfg, warpsources)
        catch
            try
                ft_topoplotER(cfg, warpsources)
            catch
                ft_topoplotTFR(cfg, warpsources)
            end
        end
        
        colorbar;
        lim=max(abs(warpsources.topo(:,currsrc)));
        lim=lim+(lim/100*10);
        caxis([-lim lim]); colorbar('off');
        colormap(bt_colorscheme('topography'));freezeColors;
        
    else % Notify user layout is preferred, if relevant
        msg2 = ['Warning: As no layout was specified in cfg.layout,',newline,...
            'the sources'' topography will not be plotted.',newline,...
            'For MEG and EEG it is recommended to specify',newline,...
            'a layout so that the warping source can be selected',newline...
            'based on activity in regions of interest.'];
        text(0,0.5,msg2,'FontName','Arial','FontSize',14,'Color',[1 0.45 0]); axis off
        if src_ind == 1
            warning('As no layout was specified in cfg.layout, the sources'' topography will not be plotted. For MEG and EEG it is recommended to specify a layout so that the warping source can be selected based on activity in regions of interest.');
        end
    end
    
    title({['Warping source ', num2str(currsrc),' (',num2str(src_ind) '/' num2str(size(srcrank,1)),')']},'FontName','Arial','FontSize',14)
    
    % Time-frequency plot
    subplot(8,2,[2 4 6 8]);
    pcolor(fspecinfo.time(mintime_ind:maxtime_ind),fspecinfo.freq,powtf(:,:,src_ind));
    shading interp
    ylim([minfft maxfft]);
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    warps_line = line([mintime,maxtime],[maxfreq, maxfreq],'Color',bt_colorscheme('warpingsignal'),'LineWidth',4);
    minfoi_l = line([mintime,maxtime],[minfoi, minfoi],'Color',bt_colorscheme('foi_borders'),'LineWidth',1);
    maxfoi_l = line([mintime,maxtime],[maxfoi, maxfoi],'Color',bt_colorscheme('foi_borders'),'LineWidth',1);
    text((maxtime+mintime)/2,maxfreq+1,[num2str(maxfreq),' Hz'],'FontSize',14,'fontweight', 'bold','Color','k')
    legend(warps_line','Warping signal')
    set(gca,'FontSize',14);
    title(sprintf('Warping source %d | Warping signal: %0.3f power at %0.2fHz',currsrc,(srcrank(src_ind,4))/maxp,maxfreq),'FontSize',14)
    colormap(bt_colorscheme('ephys_tfr'));freezeColors;
    
    % Dominant oscillation plot
    subplot(8,2,[12 14 16]);
    
    % Mark max pow and frequency of interest window
    xvec = fspecinfo.freq;
    yvec = squeeze(pspec(:,src_ind,1));
    maxpowloc = nearest(xvec,maxfreq); %Identify the warping signal location (max power)
    plot(xvec,yvec./maxp,'LineWidth',5,'Color',bt_colorscheme('foi_borders'));hold on
    w_mrk = plot(xvec(maxpowloc),yvec(maxpowloc)./maxp,'o','MarkerSize',15,'LineWidth',4,'Color',bt_colorscheme('warpingsignal'));
    w_ln = plot([xvec(maxpowloc) xvec(maxpowloc)],[-0.03,yvec(maxpowloc)/maxp],'LineWidth',3,'Color',bt_colorscheme('warpingsignal'));
    text(xvec(maxpowloc)-(xvec(2)-xvec(1))/2,1.075,[num2str(maxfreq),' Hz'],'Color','k','fontweight', 'bold','FontSize',14);
    
    try % I got you, old Matlab version users
        xline(minfoi,'Color',bt_colorscheme('foi_borders'),'LineWidth',1.5);
        xline(maxfoi,'Color',bt_colorscheme('foi_borders'),'LineWidth',1.5);
    catch
        hline(minfoi);
        hline(maxfoi);
    end
    
    set(gca,'FontSize',14);
    legend(w_mrk,'Warping signal',['Warping signal (',num2str(maxfreq),' Hz'])
    xlim([minfft maxfft]);
    ylim([-0.03 1.03])
    xlabel('Frequency (Hz)')
    ylabel('Power (a.u.)'); hold off
       
    % Adapt figure based on keypress
    keydown = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    [~,~,keyCode] = KbCheck;
    key = KbName(find(keyCode));
    
    if exist('inst','var')
        delete(inst)
    end
    
    if (keydown == 0) %grab the source if click
        src_oi=src_ind;
        if exist('ann','var')
            delete(ann)
        end
        ann = annotation('textbox',[0.18 0.19 .5 .5],'String',[' Warping source ', num2str(currsrc) ,' is highlighted'],'FitBoxToText','on','Color',bt_colorscheme('warpingsignal'),'FontSize',20,'fontweight','bold','LineWidth',1.5,'EdgeColor',[0.1 0.1 0.1]);axis off
    elseif value == 28 %go back previous source if press back arrow
        src_ind = src_ind-1;
    elseif value == 29 %go back previous source if press back arrow
        src_ind = src_ind+1;
        if src_ind > numel(srcrank) % make sure source cannot go out of bounds
            src_ind = 1;
        end
    elseif sum(strcmp(key,'q'))>=1 || sum(strcmp(key,'space'))>=1 || sum(strcmp(key,'return'))>=1 %stop the loop if it is not necessesary to keep visualising
        fprintf('Warping signal will be the %0.2fHz phase in warping source %d.',maxfreq,currsrc)
        src_ind = (numel(srcrank))+1;
        close(f1)
        finish = 1;
    end
end



%% Save basic info
bt_source{1} = src_oi;                             % Warping source which contains the warping signal
bt_source{2} = phs(src_oi);                        % Phase of all frequencies in this warping source
bt_source{3} = warpsources;                        % Warping sources data
bt_source{4} = srcrank(src_oi,:);                  % Time freq data of chosen warping source
bt_source{5} = fspecinfo;                          % FFT time and frequency vector
bt_source{6} = cutmethod;                          % Applied cutting method