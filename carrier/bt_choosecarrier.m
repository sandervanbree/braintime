function [bt_carrier] = bt_choosecarrier(config, fft_chans, channels)
% Display the top channels sorted by the characteristics of interest
% (done by bt_analyzechannels). After carrier selection, the phase of
% the oscillation with the most power in the frequency range of interest
% will be extracted. The toolbox designates this phase vector as
% representative of brain time.
%
% Use:
% [bt_carrier] = bt_choosecarrier(cfg,fft_channels,channels)
%
% Input Arguments:
% config
%   - layout         % A layout file to plot channel data at the sensor
%                    % level.
%                    %
% fft_channels       % Data structure with time frequency characteristics
%                    % of all channels as obtained by bt_analyzechannels.
%                    %
% channels           % FieldTrip channel data structure that contains
%                    % ICA components, virtual channels, or LFP time
%                    % series. One channel contains the carrier that will
%                    % be designated as brain time in this function.
%                    %
% Output:            %
% bt_carrier         % Data structure with: chosen carrier, its
%                    % time frequency information, and config details

%% Get basic info
chanrank = fft_chans{1};                             % Time freq data of top channels
mintime_ind = fft_chans{2}(1);                       % Index of start time of interest (differs for cutartefact)
maxtime_ind = fft_chans{2}(2);                       % Index of end time of interest
minfft = fft_chans{3}(1);                            % Lowest analyzed freq in FFT
maxfft = fft_chans{3}(2);                            % Highest analyzed freq in FFT
minfoi = fft_chans{4}(1);                            % Lowest possible warping frequency
maxfoi = fft_chans{4}(2);                            % Highest possible warping frequency
fspecinfo = fft_chans{5};                            % FFT time and frequency vector
powtf = fft_chans{6};                                % Power spectrum of channels
pspec = fft_chans{7};                                % Power spectrum averaged across trials
phs = fft_chans{8};                                  % Phase of all channels for all trials
cutmethod = fft_chans{9};                            % Applied cutting method

%% Plot top channels
% Find out sorting type
if isfield(config,'layout')
    sorttype = 'topography of interest';
else
    sorttype = 'maximum power at frequency of interest';
end

% Prepare instruction message
msg = ['The channels are sorted based on ',sorttype,newline,...
    'Please pick the channel containing your carrier, based on:',newline,...
    '(1) time frequency characteristics',newline,...
    '(2) topography of interest (if relevant)',newline,...
    '(3) variance explained (for ICA components, lower channel value means higher r^2)',newline,newline,...
     'INSTRUCTIONS:',newline,...
     'Press \leftarrow or \rightarrow arrow to see the previous or next channel',newline,...
     'Click on a channel to highlight it, then keep browsing or press ''Q'' to select highlighted channel'];

% Prepare figure
channeloi = -2; %channel of interest
caxislim=max(max(max(powtf(:,:,:)))); %establish axis limit
finish = 0;
chanind=1;
f1 = figure('units','normalized','outerposition',[0 0 1 1]);hold on;

% Create loop to browse through channels 
while finish==0
    if chanind < 1 % make sure channel cannot go out of bounds
        chanind = 1;
    elseif chanind>size(chanrank,1)
        chanind = size(chanrank,1);
    end
    
    % Establish current channel and its maximum frequency
    currchan = chanrank(chanind);
    maxfreq = chanrank(chanind,2);
  
   
    subplot(8,2,[1 3]);
    text(0,0.5,msg,'FontName','Arial','FontSize',14); axis off
    
    subplot(8,2,[7 9 11 13 15]);
    % Only plot topography if layout is specified
    if isfield(config,'layout')
        % channel topography
        cfg           = [];
        cfg.component = currchan; % specify the channel(s) that should be plotted
        cfg.layout    = config.layout; % specify the layout file that should be used for plotting
        cfg.comment   = 'no';
        
        % "Channels" could be either ICA comps or another data structure
        try
            ft_topoplotIC(cfg, channels)
        catch
            try
                ft_topoplotER(cfg, channels)
            catch
                ft_topoplotTFR(cfg, channels)
            end
        end
        
        colorbar;
        lim=max(abs(channels.topo(:,currchan)));
        lim=lim+(lim/100*10);
        caxis([-lim lim])
    else
        delete(f1);
        if chanind == 1
            warning('As no layout was specified in cfg.layout, the channels''s topography will not be plotted. For MEG and EEG it is recommended to specify a layout so that the brain time channel can be chosen based on activity in regions of interest.');
        end
    end
    
    title({['Channel ', num2str(currchan),' (',num2str(chanind) '/' num2str(size(chanrank,1)),')']},'FontName','Arial','FontSize',14)
    
    % Time-frequency plot
    subplot(8,2,[2 4 6 8]);
    pcolor(fspecinfo.time(mintime_ind:maxtime_ind),fspecinfo.freq,powtf(:,:,chanind));
    shading interp
    ylim([minfft maxfft]);
    caxis([0 caxislim])
    xlabel('Time (s)')
    ylabel('Frequency')
    yl = yline(maxfreq,':',[num2str(maxfreq) 'Hz'],'LineWidth',4,'Color','r');
    rec1 = rectangle('Position',[fspecinfo.time(mintime_ind),minfoi,(fspecinfo.time(maxtime_ind)-fspecinfo.time(mintime_ind)),maxfoi-minfoi],'EdgeColor','k');
    legend(yl','Carrier oscillation')
    set(gca,'FontSize',14);
    title(sprintf('Channel %d | Carrier: %0.3f power at %0.2fHz',currchan,chanrank(chanind,4),maxfreq),'FontSize',14)

    % Dominant oscillation plot
    subplot(8,2,[12 14 16]);
    
    % Mark max pow and frequency of interest window
    xvec = fspecinfo.freq;
    yvec = squeeze(pspec(:,chanind,1));
    ymax = max(max(squeeze(pspec(:,:,1))));
    maxpowloc = nearest(xvec,maxfreq); %Identify the carrier location (max power) 
    plot(xvec,yvec,'LineWidth',5);hold on
    mrk = plot(xvec(maxpowloc),yvec(maxpowloc),'o','MarkerSize',15,'LineWidth',4,'Color','r');
    mx = plot([xvec(maxpowloc) xvec(maxpowloc)],[min(yvec),ymax],'LineWidth',3,'Color','r');
    text(xvec(maxpowloc)+(xvec(2)-xvec(1))/5,(ymax+min(yvec))/2,[num2str(maxfreq),' Hz'],'Color','red','FontSize',14);
    rec2 = rectangle('Position',[minfoi,min(yvec),maxfoi-minfoi,ymax-min(yvec)],'FaceColor',[0, 0, 0, 0.15]);
    set(gca,'FontSize',14);
    legend(mrk,'Carrier oscillation',['Carrier oscillation (',num2str(maxfreq),' Hz'])
    xlim([minfft maxfft]);
    ylim([0 caxislim])
    xlabel('Frequency')
    ylabel('Power'); hold off
    
    % Adapt figure based on keypress
    keydown = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    [~,~,keyCode] = KbCheck;
    key = KbName(find(keyCode));

    if (keydown == 0) %grab the channel if click
        channeloi=chanind;
        if exist('ann','var')
            delete(ann)
        end
        ann = annotation('textbox',[0.18 0.19 .5 .5],'String',['Warping source ', num2str(currchan) ,' is highlighted'],'FitBoxToText','on','Color','red','FontSize',20);axis off
    elseif value == 28 %go back previous channel if press back arrow
        chanind = chanind-1;
    elseif value == 29 %go back previous channel if press back arrow
        chanind = chanind+1;
        if chanind > numel(chanrank) % make sure channel cannot go out of bounds
            chanind = 1;
        end
    elseif strcmp(key,'q') %stop the loop if it is not necessesary to keep visualising
        fprintf('Carrier will be the %0.2fHz phase in channel %d.',maxfreq,currchan)
        chanind = (numel(chanrank))+1;
        close(f1)
        finish = 1;   
    end
end

%% Save basic info
bt_carrier{1} = channeloi;                             % Channel which contains the carrier
bt_carrier{2} = phs(channeloi);                        % Phase of all frequencies in this channel
bt_carrier{3} = channels;                              % Channel data
bt_carrier{4} = chanrank(channeloi,:);                 % Time freq data of chosen channel
bt_carrier{5} = fspecinfo;                             % FFT time and frequency vector
bt_carrier{6} = cutmethod;                             % Applied cutting method
bt_carrier{7} = [mintime_ind, maxtime_ind];            % Index of start and end time of interest (differs for cutartefact) 
bt_carrier{8} = 'bt_choosecarrier';                    % Save carrier choosing method
