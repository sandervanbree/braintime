function [warpsigGED] = ged_foiComp(data,warpfreq)
% This function extracts a single timeseries at around the warping
% frequency based on all the warping sources using General
% Eigendecomposition (GED). This yields one timeseries per trial that is
% a holistic estimation of the warping signal, as its estimation involves
% all warping sources.
%
% Use:
% [warpfreqGED] = ged_foiComp(data,wfreq)
%
% Input Arguments:
%
% data               % Warping sources structured in a timelocked
%                    % FieldTrip format
%                    %
% wfreq              % The selected warping frequency
%                    %
% Output:            %
% warpsigGED         % 2D matrix of GED data (trial x time)
%
% These functions were written by Mike X Cohen (see Cohen 2017, eLife;
% doi: 10.7554/eLife.21792), and modified by Luca Kolibius and
% Sander van Bree.

%% Parameters
settings.trials  = 1:size(data.trial,1);                     % number of trials
settings.foi     = warpfreq;                                 % warping frequency
settings.fs      = round(1/(data.time(2)-data.time(1)));     % sampling rate
settings.visopt  = 1;                                        % visualize results? Warning: only 3 example trials will be plotted

% Full width at half maximum around the warping frequency (higher fwhm will include a more frequencies)
% Use exponential function to determine FWHM
k = 0.15;                                % Exponent
Tc = 0.5;                                % Horizontal asymptote
T0 = 3.5;                                % Vertical asymptote
fwhm = Tc+(T0-Tc).*exp(-k.*warpfreq);    % Scale FWHM with frequency
settings.fwhm   = fwhm;                                        

GEDinput = data.trial;                                       % Rename
trlLen   = size(GEDinput,3);                                 % trial length (number of samples)

warpsigGED = zeros(size(data.trial,1),trlLen);               % pre-allocate data

for trl = settings.trials
    %% GED to identify warping frequency component
    % Warp frequency covariance
    currTrl = squeeze(GEDinput(trl,:,:));
    trlLFPfilt = filterFGx(currTrl,settings.fs,settings.foi,settings.fwhm); % channel x time filtered around warping frequency
    trlLFPfilt = bsxfun(@minus,trlLFPfilt,mean(trlLFPfilt,2));              % demeaning
    wfreqcov  = (trlLFPfilt*trlLFPfilt')/trlLen;                            % calculate covariance matrix (normalized by length)
    
    % broadband covariance
    trlLFP = squeeze(GEDinput(trl,:,:));                                    % same for broadband signal
    bbcov = (trlLFP*trlLFP') / trlLen;
    
    % GED
    [evecsT,evals] = eig(wfreqcov,bbcov);   % compute generalized eigendecomposition. this is the vector that maximizes similarity of wfreqcov to bbcov
    
    % find best component and compute filter projection
    [~,maxcomp] = sort(diag(evals));        % sort eigenvalues so the last one is the highest
    wfreqmap    = wfreqcov*evecsT;
    wfreqmap    = wfreqmap(:,maxcomp(end));
    
    % fix sign of map (max is positive)
    [~,maxe] = max(abs(wfreqmap));
    wfreqmap = wfreqmap * sign(wfreqmap(maxe));
    
    % wfreq time series component
    wfreqcomp = trlLFPfilt' * evecsT(:,maxcomp(end)); % weigh the filtered data with the highest eigenvector
    
    % fix sign of time series according to sign of correlation with ephys data
    wfreqcomp = wfreqcomp * sign(corr(wfreqcomp,filterFGx(trlLFP(maxe,:),settings.fs,settings.foi,9)'));
    
    % Output matrix
    warpsigGED(trl,:) = wfreqcomp;
    
    %% Optional visualization
    if settings.visopt == 1
        if trl <= 3 % Only plot 3 example trials
            figure;
            subplot(311)
            plot(data.time,wfreqcomp, 'linew', 2, 'color', [0 0 0 ])
            xticks('')
            ylabel(' \muV (?)')
            xlim([data.time(1) data.time(end)])
            title('Warping frequency Component')
            mAx = gca;
            mAx.FontWeight = 'bold';
            mAx.FontSize   = 20;
            
            subplot(312);
            meanFoi = mean(trlLFPfilt,1);
            plot(data.time,meanFoi, 'linew', 2, 'color', [0 0 0]);
            xticks('')
            ylabel('amp')
            xlim([data.time(1) data.time(end)])
            title('Average warping frequency over channels')
            mAx = gca;
            mAx.FontWeight = 'bold';
            mAx.FontSize   = 20;
            
            subplot(313); hold on
            title('All channels')
            for chan = 1:size(trlLFP,1)
                plot(data.time,(trlLFPfilt(chan,:)+0.01*chan), 'color', [0 0 0 0.5])
            end
            ylim([0 0.17])
            xlim([data.time(1) data.time(end)])
            xlabel('Time (s)')
            yticks('')
            ylabel({'Channel', '(stacked)'})
            mAx = gca;
            mAx.FontWeight = 'bold';
            mAx.FontSize   = 20;
        end
    end
end
end
