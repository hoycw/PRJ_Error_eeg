function norm_filt = fn_bsln_ft_filtered(filt, bsln_lim, bsln_type, n_boots)
%% Baseline correct one TFR based on bsln_lim epoch, both from ft_freqanalysis
% INPUTS:
%   filt [ft dataset] - full output of ft_preprocessing (no .powspctrm!)
%   bsln_lim [int, int]- 2 int array of TIME indices for [start, end] of baseline period
%   bsln_type [str]   - type of baseline to implement
%       'zboot'  = pool all baselines, bootstrap means+SDs, z-score all
%           trials to the mean+SD of the bootstrapped distribution
%       'zscore' = subtract mean and divide by SD
%       'demean' = subtract mean
%       'my_relchange' = subtract mean, divide by mean (results in % change)
%   n_boots [int] - number of bootstrap iterations
% OUTPUTS:
%   bslnd_tfr [ft dataset] - same tfr but baseline corrected

[~, app_dir] = fn_get_root_dir();
addpath([app_dir 'fieldtrip/']);
ft_defaults
rng('shuffle'); % seed randi with time

% Select baseline data
cfgs = [];
cfgs.latency = bsln_lim;
bsln_filt = ft_selectdata(cfgs,filt);

norm_filt = filt;
for ch = 1:numel(filt.label)
    % Select this channel
    fprintf('\tBaseline correcting %s (%i / %i)\n',filt.label{ch},ch,numel(filt.label));
    cfgs = [];
    cfgs.channel = filt.label(ch);
    ch_bsln = ft_selectdata(cfgs,bsln_filt);
    ch_filt = ft_selectdata(cfgs,filt);
    
    % Create bootstrap distribution if necessary
    if strcmp(bsln_type,'zboot')
        sample_means = NaN([1 n_boots]);
        sample_stds  = NaN([1 n_boots]);
        for boot_ix = 1:n_boots
            % Select a random set of trials (sampling WITH REPLACEMENT!)
            shuffle_ix = randi(numel(filt.trial),[1 numel(filt.trial)]);
            % Select baseline data and compute stats
            bsln_data = horzcat(ch_bsln.trial{shuffle_ix});
            sample_means(boot_ix) = nanmean(bsln_data(:));
            sample_stds(boot_ix)  = nanstd(bsln_data(:));
        end
    end
    
    % Perform Baseline Correction
    for t = 1:numel(filt.trial)
        trial    = ch_filt.trial{t};
        trl_bsln = ch_bsln.trial{t};
        switch bsln_type
            case 'zboot'
                norm_filt.trial{t}(ch,:) = (trial-mean(sample_means))/mean(sample_stds);
            case 'zscore'
                norm_filt.trial{t}(ch,:) = (trial-nanmean(trl_bsln))/nanstd(trl_bsln);
            case 'demean'
                norm_filt.trial{t}(ch,:) = trial-nanmean(trl_bsln);
            case 'my_relchange'
                norm_filt.trial{t}(ch,:) = (trial-nanmean(trl_bsln))/nanmean(trl_bsln);
            otherwise
                error(['Unknown bsln_type: ' bsln_type])
        end
    end
end

end
