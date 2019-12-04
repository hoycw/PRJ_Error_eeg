function SBJ04c_ERP_grp_stats_GLM(SBJs,proc_id,an_id,stat_id,varargin)
% Compute grand average group ERP from SBJ ERPs:
%   Re-align data to event, select channels and epoch, filter, average, run stats, save
% INPUTS:
%   SBJs [cell array] - ID list of subjects to run
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:
%   grp_erp [cell] - cell array with outputs of ft_timelockgrandaverage for each condition
%   NOT stat [ft struct] - output of ft_timelockstatistics, not done yet

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
%         elseif strcmp(varargin{v},'write_report')
%             write_report = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ischar(save_fig); save_fig = str2num(save_fig); end
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' stat_id '/' an_id '/'];
if ~exist(fig_dir,'dir') && save_fig
    mkdir(fig_dir);
end

%% Load Data 
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.measure,'p2p')
    error('peak latencies need jackknife and p2p!');
end

% Select Conditions of Interest
[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles(conditions);

% Load example data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{1} '_vars.m'];
eval(SBJ_vars_cmd);
tmp = load([root_dir 'PRJ_Error_eeg/data/',SBJs{1},'/04_proc/',SBJs{1},'_',an_id,'.mat'],'roi');
ch_list  = tmp.roi.label;
cfg = []; cfg.latency = st.stat_lim;
st_roi = ft_selectdata(cfg, tmp.roi);
time_vec = st_roi.time{1};

% Load Data
erps = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(time_vec)]);
for s = 1:length(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    tmp = load([SBJ_vars.dirs.proc,SBJs{s},'_',an_id,'.mat']);
    load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);
    bhvs{s} = tmp.bhv; rois{s} = tmp.roi; w2s{s} = tmp.w2;
    clear SBJ_vars tmp
end

%% Compute Grand Average Group ERPs
cfg_trim = []; cfg_trim.latency = st.stat_lim;
st_data = nan([numel(SBJs) numel(w2{1}.label) numel(w2{1}.time)]);
for s = 1:numel(SBJs)
    tmp = ft_selectdata(cfg, rois{s});
    st_data(s,:,:) = mean(w2{s}.zscore;
end
if any(isnan(st_data(:))); error('NaN in ANOVA data!'); end

% grp_erp = cell(size(cond_lab));
% for cond_ix = 1:numel(cond_lab)
%     grp_erp{cond_ix} = ft_timelockgrandaverage(cfg_gavg,roi_erps{:,cond_ix});
% end

%% Run Statistics

%% Save Results
data_out_fname = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',conditions,'_',an_id,'.mat');
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'-v7.3','grp_erp','SBJs');

end
