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
FFT_phs = fft_sources{7};                              % Phase of all warping sources for all trials
cutmethod = fft_sources{8};                            % Applied cutting method
rankmethod = fft_sources{9};                           % Applied ranking method
asymmidx = fft_sources{10};                            % Asymmetry index per warping source
asymmidx_t = fft_sources{11};                          % Asymmetry index t-stat per warping source
wavshap = fft_sources{12};                             % Waveshape per warping source

mintime_ind = nearest(fspecinfo.time,mintime);         % Index of start time of interest (differs for cutartefact)
maxtime_ind = nearest(fspecinfo.time,maxtime);         % Index of end time of interest
    
%% Plot ranked warping sources
% Find out ranking type
if strcmp(rankmethod,'maxpow')
    ranktype = 'maximum power.';
elseif strcmp(rankmethod,'templatetopo')
    ranktype = 'maximum power and correlation to template topography.';
end

% Prepare instruction message
msg = ['The warping sources are ranked by ',newline,...
    ranktype,newline,newline,...
    'Select a warping source based on its:',newline,...
    '(1) time frequency characteristics',newline,...
    '(2) topography (if applicable)',newline,...
    '(3) variance explained (lower ICA',newline...
    'components explain more variance)',newline,newline,...
    'Browse:       \leftarrow or \rightarrow arrow',newline,...
    'Highlight:     Click component',newline,...
    'Quit:            Q, X, or space'];

% Prepare figure
src_oi = -2; %source of interest
maxp=max(max(max(powtf(:,:,:)))); %establish maximum power. This will be used to limit axes and normalize power to
% Sanity check
if isnan(maxp)
    error(['The time frequency analysis of the warping source analysis has yielded only NaNs.'...
        ' Please check whether (1) the correct time window was used for the analysis, and (2)'...
        ' whether your data contains sufficient additional time before and after your time window.']);
end
finish = 0;
src_ind=1;
f1 = figure;hold on;

% Adapt figure
bt_figure;

if isfield(config,'quickselect') && config.quickselect == 1
    finish = 1;
    src_oi=src_ind;
    close(f1);
end

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
    
    subplot(11,11,[1 35]); % Plot text
    inst = text(0,0.5,msg,'FontName','Arial','FontSize',15); axis off
    
    % Time-frequency plot
    subplot(11,11,[67 115])
    pcolor(fspecinfo.time(mintime_ind:maxtime_ind),fspecinfo.freq,powtf(:,:,src_ind));
    shading interp
    ylim([minfft maxfft]);
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    warps_line = line([mintime,maxtime],[maxfreq, maxfreq],'Color',bt_colorscheme('warpingsignal'),'LineWidth',4);
    minfoi_l = line([mintime,maxtime],[minfoi, minfoi],'Color',bt_colorscheme('foi_borders'),'LineWidth',1);
    maxfoi_l = line([mintime,maxtime],[maxfoi, maxfoi],'Color',bt_colorscheme('foi_borders'),'LineWidth',1);
    text((maxtime+mintime)/2,maxfreq+1,[num2str(maxfreq),' Hz'],'FontSize',14,'fontweight', 'bold','Color','k')
    l1 = legend(warps_line','Warping signal'); 
    l1.FontSize = 14;
    set(gca,'FontSize',14); 
    title(sprintf('Warping source %d | Warping signal: %0.3f power at %0.2fHz',currsrc,(srcrank(src_ind,4))/maxp,maxfreq),'FontSize',14)
    colormap(bt_colorscheme('ephys_tfr'));freezeColors;
    
    % Topography
    sp = subplot(11,11,[3 39]);
    % Only plot topography if layout is specified
    if isfield(config,'layout')
        % warping source topography
        cfg           = [];
        cfg.figure    = sp;
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
    
    % Dominant oscillation plot
    subplot(11,11,[73 121]);
    
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
    l2 = legend(w_mrk,'Warping signal',['Warping signal (',num2str(maxfreq),' Hz']);
    l2.FontSize = 14;
    xlim([minfft maxfft]);
    ylim([-0.03 1.03])
    xlabel('Frequency (Hz)')
    ylabel('Power (a.u.)'); hold off

    % Plot asymmetry across cycles
    asymm_prc = prctile(abs(asymmidx(:)),95);
    asymm_lim = [-asymm_prc asymm_prc];
    
    try  % This code only makes sense if waveshape and asymmetry have been calculated
        sp2 = subplot(11,11,[7 53]);hold on; % Plot asymmetry across all cycles
        bt_rainplot(asymmidx(src_ind,:),[],16,30,asymm_lim);
        tval = mean(asymmidx_t(src_ind,:));
        title(['Asymmetry: ',num2str(round(tval,3)),' (t)']);
        
        sp3 = subplot(11,11,[10,55]);hold on; % Plot waveshape
        bt_wavplot(wavshap(src_ind,:),2);
        title('Average waveshape');
    catch
    end
    
    % Adapt figure based on keypress
    keydown = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    
    % Delete plots
    if exist('inst','var')
        delete(inst)
    end   
    if exist('sp2','var')
        delete(sp2)
    end
    if exist('sp3','var')
        delete(sp3)
    end
    
    if (keydown == 0) %grab the source if click
        src_oi=src_ind;
        if exist('ann','var')
            delete(ann)
        end
        ann = annotation('textbox',[0.18 0.075 .45 .5],'String',['Warping source ', num2str(currsrc) ,' is highlighted'],'FitBoxToText','on','Color',bt_colorscheme('warpingsignal'),'FontSize',20,'fontweight','bold','LineWidth',1.5,'EdgeColor',[0.1 0.1 0.1]);axis off
    elseif value == 28 %go back previous source if press back arrow
        src_ind = src_ind-1;
    elseif value == 29 %go back previous source if press back arrow
        src_ind = src_ind+1;
        if src_ind > numel(srcrank) % make sure source cannot go out of bounds
            src_ind = 1;
        end
    elseif value == 113 || value == 87  || value == 120 || value == 88 || value == 32 %stop the loop if it is not necessesary to keep visualising
        fprintf('Warping signal will be the %0.2fHz phase in warping source %d.',maxfreq,currsrc)
        src_ind = (numel(srcrank))+1;
        close(f1)
        finish = 1;
    end
end

%% Get Generalized Eigendecomposition (GED) phase for warping frequency
% A timelock format is well-geared for Mike X Cohen's functions.
cfg                     = [];
cfg.keeptrials          = 'yes';
cfg.removemean          = 'no';
GEDdata                 = ft_timelockanalysis(cfg, warpsources);
warpfreq                = srcrank(src_oi,2);

% Call Mike X Cohen's GED functions
[warpsigGED] = ged_foiComp(GEDdata,warpfreq);

% Find time bins of the FFT's start and end time
[~,mintime_fft_ind] = min(abs(fspecinfo.time(1)-warpsources.time{1}));
[~,maxtime_fft_ind] = min(abs(fspecinfo.time(end)-warpsources.time{1}));

% Crop to time window of interest
warpsigGED = warpsigGED(:,mintime_fft_ind:maxtime_fft_ind);

% Pre-allocate
GED_phs = zeros(size(warpsigGED,1),size(warpsigGED,2));

% Extract phase of Hilbert transformed GED signal
for trl = 1:size(warpsigGED,1)
    GED_phs(trl,:) = angle(hilbert(warpsigGED(trl,:)));
end

%% Save basic info
bt_source{1} = src_oi;                             % Warping source which contains the warping signal
%               ^ this is the source index, not channel name!
bt_source{2} = FFT_phs(src_oi);                    % Phase of all frequencies in this warping source
bt_source{3} = warpsources;                        % Warping sources data
bt_source{4} = srcrank(src_oi,:);                  % Time freq data of chosen warping source
bt_source{5} = fspecinfo;                          % FFT time and frequency vector
bt_source{6} = cutmethod;                          % Applied cutting method
bt_source{7} = wavshap(src_oi,:);                  % Waveshape of warping signal
bt_source{8} = mean(asymmidx(src_oi,:));           % Asymmetry index of warping signal
bt_source{9} = GED_phs;                            % Phase of warping signal as estimated using GED