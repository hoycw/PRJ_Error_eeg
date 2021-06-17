%% Plot single electrode topography
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

addpath([app_dir 'fieldtrip/']);
ft_defaults


elec_lab  = 'Fz';
tmp_SBJ   = 'EEG01';
tmp_an_id = 'ERP_all_F2t1_dm2t0_fl05t20';

load([root_dir 'PRJ_Error_eeg/data/' tmp_SBJ '/04_proc/' tmp_SBJ '_' tmp_an_id '.mat']);

erp = ft_timelockanalysis([],roi);
e_ix = find(strcmp(erp.label,elec_lab));
erp.avg(setdiff(1:numel(erp.label),e_ix),:) = nan;
erp.avg(e_ix,:) = 10;

cfgp = [];
cfgp.layout           = 'biosemi64.lay';
cfgp.comment          = 'no';
cfgp.marker           = 'off';
cfgp.highlight        = 'on';
cfgp.highlightsymbol  = '*';
cfgp.highlightsize    = 10; % default = 6
cfgp.highlightcolor   = [1 0 0];
cfgp.highlightchannel = {elec_lab};
cfgp.style            = 'blank';
ft_topoplotER(cfgp, erp);
% title(['Electrode ' elec_lab]);
