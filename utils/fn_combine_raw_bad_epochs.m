function bad_epochs = fn_combine_raw_bad_epochs(SBJ) 
%% Combine bad epoch indices after concatenating across blocks
% INPUTS:
%   SBJ [str] - name of SBJ
% OUTPUTS:
%   bad_epochs [Nx2 array] - [start stop] indices of N bad epochs for concatenated dataset

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

% Load bad epochs from SBJ00 viewing
if numel(SBJ_vars.block_name)>2; error('Sheila code does not work for more than 2 blocks!'); end
if numel(SBJ_vars.block_name)==1
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

% Correct time stamps for concatenation
for b_ix = 2:numel(SBJ_vars.block_name)
    initial_bad_epochs{b_ix} = initial_bad_epochs{b_ix} + SBJ_vars.endsample{b_ix-1}/SBJ_vars.origsample*SBJ_vars.resample;  
end

% Combine bad epochs
bad_epochs = vertcat(initial_bad_epochs{1},initial_bad_epochs{2});
% functional but poor logic Sheila code:
    % for b_ix = 2:numel(SBJ_vars.block_name)
    %     initial_bad_epochs{b_ix} = vertcat(initial_bad_epochs{b_ix-1},initial_bad_epochs{b_ix});
    %     bad_epochs = initial_bad_epochs{b_ix};
    % end

end
