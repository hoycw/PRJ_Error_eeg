function bad_epochs = fn_combine_raw_bad_epochs(SBJ) 

if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Processing variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
if numel(SBJ_vars.block_name)==1;
    block_suffix = '';
    in_fname = [SBJ_vars.dirs.events SBJ '_raw_bad_epochs' block_suffix '.mat'];
    load(in_fname);
    initial_bad_epochs = bad_epochs;
else
    for block = 1:numel(SBJ_vars.block_name)
        block_suffix = ['_' SBJ_vars.block_name{block}];
        in_fname = [SBJ_vars.dirs.events SBJ '_raw_bad_epochs' block_suffix '.mat'];
        load(in_fname)
        initial_bad_epochs{block} = bad_epochs;
    end
end
for b_ix = 2:numel(SBJ_vars.block_name)
    initial_bad_epochs{b_ix} = initial_bad_epochs{b_ix} + SBJ_vars.endsample{b_ix-1}/origsample*resample;  
end
for b_ix = 2:numel(SBJ_vars.block_name)
    initial_bad_epochs{b_ix} = vertcat(initial_bad_epochs{b_ix-1},initial_bad_epochs{b_ix});
    bad_epochs = initial_bad_epochs{b_ix};
end

