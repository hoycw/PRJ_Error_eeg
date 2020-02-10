%% Create cfg.trl
cfg = [];
cfg.dataset = SBJ_vars.dirs.raw_filename;
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.eventtype = 'gui';
cfg = ft_definetrial(cfg)

%% Build trigger channel
trigs = zeros([1 size(data.trial{1},2)]);
for ev_ix = 1:size(cfg.trl,1);
    trigs(cfg.trl(ev_ix,1):cfg.trl(ev_ix,2)) = cfg.trl(ev_ix,4);
end

%% Plot the triggers
plot(data.time{1},trigs);

%% Extract stim onsets to get ITI
% Get all ITIs
iti = zeros([size(cfg.trl,1) 1]);
last_stim = 0;
for ev_ix = 2:size(cfg.trl,1)
    if cfg.trl(ev_ix,4)==1
        iti(ev_ix) = cfg.trl(ev_ix,1)-last_stim;
        last_stim = cfg.trl(ev_ix,1);
    end
end
% Remove non-stim ITIs and between block ITIs
iti(iti==0) = [];
bw_block_times = iti(iti>5*data.fsample)/data.fsample;
iti(iti>5*data.fsample) = [];
iti_s = iti/data.fsample;
