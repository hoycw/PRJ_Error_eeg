function BHV01_plot_RT_hist(SBJ,proc_id,save_fig,varargin)
% Plot RT histogram for example SBJ

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');      fig_vis = 'on'; end
if ~exist('fig_ftype','var');    fig_ftype = 'png'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Load Data 
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);

load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

%% Plot Predictors
fig_name = [SBJ '_RT_histogram'];
figure('Name',fig_name,'Visible',fig_vis);%,'units','normalized','OuterPosition',[0 0 0.5 0.5]);
ax = gca; hold on;

% Plot histogram
histogram(bhv.rt);

% Plot Kernel Density Estimate
% [kde, kde_ix] = ksdensity(bhv.rt);
% plot(kde_ix,kde,'Color','k');
xlabel('Reaction Time (s)');
title([SBJ ' mean = ' num2str(mean(bhv.rt))],'Interpreter','none');
set(gca,'FontSize',16);

%% Save Figure
if save_fig
    % Create figure directory
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/RTs/histograms_MATLAB/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    % % Commented out because screen ratios aren't right automatically
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
