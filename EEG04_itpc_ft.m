%ITPC/TFA Script
function itpc(filename, SBJ)
addpath(genpath('/Users/SCS22/Desktop/Knight_Lab/Preprocessing_Work'));
load(filename,'data');
SBJ_vars_cmd = ['run ' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
cfg = [];             
cfg.method     = 'wavelet';                
cfg.width      = 4; 
cfg.output     = 'pow';	
cfg.foilim     = [2 30];                
cfg.toi        = 0:0.004:1; 

TFRwave = ft_freqanalysis(cfg, data);
cfg = [];
cfg.baseline     = [-0.1 0]; 
cfg.baselinetype = 'absolute'; 	        
%cfg.zlim         = [-3e-25 3e-25];
cfg.xlim = [-0.2 1];
cfg.showlabels   = 'yes';	
cfg.interactive = 'yes'
cfg.layout       = 'biosemi64.lay';
cfg.masknans = 'yes';
figure
ft_multiplotTFR(cfg, TFRwave)
save([filename 'TFRwave'], 'TFRwave')
savefig([data_out_filename_TFA 'multiplotTFR']);
%% ITPC
% frequency decomposition
% currently nonfunctional because the array is too big -- may be able to do
% per channel but not able to get this to work either -- all nans
cfg = [];
cfg.channel = 'Fp2';
cfg.method = 'wavelet';
cfg.toi    = 0:0.004:1;
cfg.output = 'fourier';
roi_freq = ft_freqanalysis(cfg, data);
% % ITPC computation
% itc = [];
% itc.label = roi_freq.label;
% itc.freq  = roi_freq.freq;
% itc.time  = roi_freq.time;
% itc.dimord = 'chan_freq_time';
% 
% F = roi_freq.fourierspctrm;   % copy the fourier spectrum
% N = size(F,1);                  % number of trials
% 
% itc.itpc = F./abs(F);           % divide by amplitude
% itc.itpc = sum(itc.itpc,1);     % sum angles
% itc.itpc = abs(itc.itpc)/N;     % take absolute value and normalize
% itc.itpc = squeeze(itc.itpc);
% figure
% subplot(2, 1, 1);
% imagesc(itc.time, itc.freq, squeeze(itc.itpc(1,:,:))); 
% axis xy
% title('inter-trial phase coherence');
% savefig([data_out_filename_TFA 'ITPCgraph']);
% save(data_out_filename_TFA,'itc');


