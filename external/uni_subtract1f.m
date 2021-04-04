function freq = uni_subtract1f(freq)

% This function estimates and subtracts the fractal property of a power 
% spectrum. The [freq] input must match a Fieldtrip formatted
% time-frequency data structure.
%
% By Benjamin Griffiths and Simon Hanslmayr

% extract parameters
f   = freq.freq;
t   = freq.time;
% cycle through channels
nchans = size(freq.powspctrm,2);

for ci = 1:nchans
    pow = nanmean(freq.powspctrm(:,ci,:,:),4);

    % exclude bandstop filtered frequencies from fit
    fidx = [];
    %fidx =find(f>48 & f<52 );
    f(fidx) = [];
    pow(:,fidx) = [];
    freq.freq(fidx) = [];
    freq.powspctrm(:,:,fidx,:) =[];
    


    % fit 1/f
    logf    = squeeze(log10(f));
    logpow  = squeeze(log10(pow));
%     logpow(isnan(logpow))=1; %added!
    beta    = zeros(size(logpow,1),1);            
    % cycle through each trial
    for trl = 1 : size(logpow,1)

        % get temporary power and frequency
        tmp_f   = logf;
        tmp_pow = logpow(trl,:);

        % iterate
        while true

            % get fit
            b = [ones(size(tmp_f))' tmp_f'] \ tmp_pow';
                        
            linft = tmp_f*b(2) + b(1);

            % subtract fit
            firstpass = tmp_pow - linft;

            % get error (the mean power of all frequencies less than zero)
            erange = mean(firstpass(firstpass<0));

            % find all postive values greater than the error
            eidx   = firstpass>abs(erange);

            % if no error
            if sum(eidx) == 0

                % recompute fit using all data
                linft = logf*b(2) + b(1);

                % subtract fit
                logpow(trl,:) = logpow(trl,:) - linft;

                % save output 
                beta(trl,1) = b(2);
                break

            % else if more than half of the frequencies have been removed
            elseif numel(tmp_f) <= numel(logf)/2

                % recompute fit using all data
                linft = logf*b(2) + b(1);

                % subtract fit
                logpow(trl,:) = logpow(trl,:) - linft;

                % save output 
                beta(trl,1) = b(2);
                break

             % otherwise, remove freqencies that exceed the error threshold 
            else
                tmp_pow = tmp_pow(eidx==0);
                tmp_f   = tmp_f(eidx==0);
            end
        end
        % Subtract linft from each time-point to get time-resolved spectrum
        for ti=1:length(t)
            tmp_tfpow(trl,ci,:,ti) = squeeze(freq.powspctrm(trl,ci,:,ti)) - 10.^linft';
        end
    end
end

% convert pow
freq.powspctrm = tmp_tfpow;
freq.fractal   = repmat(beta,[1 size(freq.powspctrm,2)]);
freq.freq = f;
