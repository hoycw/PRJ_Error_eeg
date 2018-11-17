
function itpc(filename, SBJ)
addpath(genpath('/Users/SCS22/Desktop/Knight_Lab/Preprocessing_Work'));
load(filename,'data');
cfg.output       = 'pow';
cfg.channel      = {'all'}
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:30;                         % analysis 2 to 30 Hz in steps of 2 Hz 
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -0.5:0.05:1.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
TFRhann = ft_freqanalysis(cfg,dataFIC);
cfg = [];
cfg.baseline     = [-0.5 -0.1]; 
cfg.baselinetype = 'absolute'; 
cfg.zlim         = [-3e-27 3e-27];	        
cfg.showlabels   = 'yes';	
cfg.interactive = 'yes'
cfg.layout       = 'biosemi64.lay';
figure 
ft_multiplotTFR(cfg, TFRhann);
savefig([data_out_filename_TFA 'multiplotTFR']);
cfg = [];
cfg.baseline     = [-0.5 -0.1];	
cfg.baselinetype = 'absolute';
cfg.xlim         = [0.9 1.3];   
cfg.zlim         = [-1.5e-27 1.5e-27];
cfg.ylim         = [15 20];
cfg.marker       = 'on';
figure 
ft_topoplotTFR(cfg, TFRhann);
savefig([data_out_filename_TFA 'topo']);
save(data_out_filename_TFA,'TFRhann');
% make a new FieldTrip-style data structure containing the ITC
% copy the descriptive fields over from the frequency decomposition
freq = TFRHann;
itc = [];
itc.label     = freq.label;
itc.freq      = freq.freq;
itc.time      = freq.time;
itc.dimord    = 'chan_freq_time';
F = freq.fourierspctrm;   % copy the Fourier spectrum
N = size(F,1);           % number of trials
% compute inter-trial phase coherence (itpc) 
itc.itpc      = F./abs(F);         % divide by amplitude  
itc.itpc      = sum(itc.itpc,1);   % sum angles
itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
figure
subplot(2, 1, 1);
imagesc(itc.time, itc.freq, squeeze(itc.itpc(1,:,:))); 
axis xy
title('inter-trial phase coherence');
