function [bt_comp] = bt_choosecomp(config, fft_comp, comp)
% Display the top components sorted by the characteristics of interest
% (done by bt_analyzecomps). After component selection, the phase of
% the oscillation with the most power in the frequency range of interest
% will be extracted. The toolbox designates this phase vector as
% representative of brain time.
%
% Use:
% [bt_comp] = bt_choosecomps(cfg,fft_comp,comp)
%
% Input Arguments:
% config
%   - layout         % A layout file to plot component data at the channel
%                    % level.
%                    %
% fft_comp           % Data structure with time frequency characteristics
%                    % of all components as obtained by bt_analyzecomps.
%                    %
% comp               % FieldTrip component data structure as obtained
%                    % by applying ft_componentanalysis on clock time data.
%                    %
% Output:            %
% bt_comp            % Data structure with: chosen component, its 
%                    % time frequency information, and config details
%                    % saved for later retrieval.

%% Get basic info
topcomps = fft_comp{1}; %What are the top components?
mintime_ind = fft_comp{2}(1);
maxtime_ind = fft_comp{2}(2);
minfft = fft_comp{3}(1); 
maxfft = fft_comp{3}(2);
fspecinfo = fft_comp{4}; %FFT characteristics
powtf = fft_comp{5};
pspec = fft_comp{6};
phs = fft_comp{7}; %Phase of all components
cutmethod = fft_comp{8};
cutinfo = fft_comp{9};

%% Plot top components
figure
compoi = [];%component of interest
caxislim=max(max(max(powtf(:,:,:)))); %establish axis limit
compind=1;

uiwait(msgbox({'The components are sorted based on power or topography (depending on sortmethod). Please pick one component, factoring in the component''s:';' ';...
    '(1) Topography.';...
    '(2) Time Frequency characteristics.';...
    '(3) optionally, the variance explained by the component (lower component number means higher r^2).'}))
uiwait(msgbox({'Instructions to browse through components:';' ';...
    'Press forward/back arrow to see the next/previous component';...
    'Once you have decided for one component, click on that component and press Q to quit visualization'}))

while compind <= numel(topcomps)
    currcomp = topcomps(compind);
    
    % component topography
    subplot(5,2,[1 3 5 7 9]);
    cfg           = [];
    cfg.component = currcomp; % specify the component(s) that should be plotted
    cfg.layout    = config.layout; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, comp)
    colorbar
    lim=max(abs(comp.topo(:,currcomp)));
    lim=lim+(lim/100*10);
    caxis([-lim lim])
    title({[num2str(compind) '/' num2str(size(topcomps,1)) ' Components'] ['Component ' num2str(currcomp)]})
   
    % time-frequency plot
    subplot(5,2,[2 4 6 8]);
    pcolor(fspecinfo.time(mintime_ind:maxtime_ind),fspecinfo.freq,powtf(:,:,currcomp));
    shading interp
    ylim([minfft maxfft]);
    caxis([0 caxislim])
    title(sprintf('Power at %0.3fHz: %0.3f',topcomps(currcomp,2),topcomps(currcomp,4)))
    
    % dominant oscillation plot
    subplot(5,2,10);
    plot(squeeze(pspec(:,currcomp,1)),fspecinfo.freq);
    ylim([minfft maxfft]);
    xlim([0 caxislim])
    
    keydown = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    [~,~,keyCode] = KbCheck;
    key = KbName(find(keyCode));
    if (keydown == 0) %grab the component if click
        compoi=currcomp;
        compind = compind-1;
    elseif value == 28 %go back previous component if press back arrow
        disp(value);
        compind = compind-2;
    elseif strcmp(key,'q') %stop the loop if it is not necessesary to keep visualising  %CHECK
      compind = (numel(topcomps))+1;
      close
    end
    compind = compind+1;
end


%% Save basic info
bt_comp{1} = compoi; %chosen component
bt_comp{2} = phs(compoi); %phase of chosen component
bt_comp{3} = comp;
bt_comp{4} = topcomps(compoi,:);
bt_comp{5} = fspecinfo;
bt_comp{6} = cutmethod;
bt_comp{7} = cutinfo;
