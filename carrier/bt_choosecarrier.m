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
%                    % saved for later retrieval.

%% Get basic info
chanrank = fft_chans{1}; %What are the top channels?
mintime_ind = fft_chans{2}(1);
maxtime_ind = fft_chans{2}(2);
minfft = fft_chans{3}(1);
maxfft = fft_chans{3}(2);
fspecinfo = fft_chans{4}; %FFT characteristics
powtf = fft_chans{5};
pspec = fft_chans{6};
phs = fft_chans{7}; %Phase of all channels
cutmethod = fft_chans{8};

%% Plot top channels
figure
channeloi = []; %channel of interest
caxislim=max(max(max(powtf(:,:,:)))); %establish axis limit
finish = 0;
chanind=1;

uiwait(msgbox({'The channels are sorted based on power or topography (depending on cfg.sortmethod). Please pick one channel, factoring in its:';' ';...
    '(1) time frequency characteristics.';...
    '(2) if relevant, topography of interest.';...
    '(3) if channels are components, the variance explained (lower channel number means higher r^2).'}))
uiwait(msgbox({'Instructions to browse through channels:';' ';...
    'Press forward/back arrow to see the next/previous channel';...
    'Once you have decided for one channel, click on that channel and press Q to quit visualization'}))

while finish==0
    if chanind < 1 % make sure channel cannot go out of bounds
        chanind = 1;
    elseif chanind>size(chanrank,1)
        chanind = size(chanrank,1);
    end
    
    currchan = chanrank(chanind);
    
    f1 = subplot(5,2,[1 3 5 7 9]);    
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
    
    title({[num2str(chanind) '/' num2str(size(chanrank,1)) ' channels']})
    
    % time-frequency plot
    subplot(5,2,[2 4 6 8]);
    pcolor(fspecinfo.time(mintime_ind:maxtime_ind),fspecinfo.freq,powtf(:,:,chanind));
    shading interp
    ylim([minfft maxfft]);
    caxis([0 caxislim])
    xlabel('Time (s)')
    ylabel('Frequency')
    title(sprintf('Channel %d | Carrier: %0.3f power at %0.2fHz',currchan,chanrank(chanind,4),chanrank(chanind,2)))
    
    % dominant oscillation plot
    subplot(5,2,10);
    plot(fspecinfo.freq,squeeze(pspec(:,chanind,1)));
    xlim([minfft maxfft]);
    ylim([0 caxislim])
    xlabel('Frequency')
    ylabel('Power')
    
    keydown = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    [~,~,keyCode] = KbCheck;
    key = KbName(find(keyCode));
    if (keydown == 0) %grab the channel if click
        channeloi=chanind;
        fprintf('Selected carrier in channel number %d. Press ''q'' to quit.',currchan)
    elseif value == 28 %go back previous channel if press back arrow
        chanind = chanind-1;
    elseif value == 29 %go back previous channel if press back arrow
        chanind = chanind+1;
        if chanind > numel(chanrank) % make sure channel cannot go out of bounds
            chanind = 1;
            disp('This is the last channel')
        end
    elseif strcmp(key,'q') %stop the loop if it is not necessesary to keep visualising
        fprintf('Carrier will be the %0.2fHz phase in channel number %d.',chanrank(chanind,2),currchan)
        chanind = (numel(chanrank))+1;
        finish = 1;
        close
    end
end

%% Save basic info
bt_carrier{1} = channeloi; %chosen carrier
bt_carrier{2} = phs(channeloi); %phase of chosen carrier
bt_carrier{3} = channels;
bt_carrier{4} = chanrank(channeloi,:);
bt_carrier{5} = fspecinfo;
bt_carrier{6} = cutmethod;
bt_carrier{7} = [mintime_ind, maxtime_ind];
