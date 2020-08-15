function [combined] = fn_concat_blocks(blocks)
%% Concatenate two fieldtrip data structures with continuous data (e.g., two task blocks)
%   Will lump any additional datasets into the first one (to keep the
%   properties of that first block and keep ft happy
% INPUTS:
%   blocks [cell array] - FT data structs to be concatenated (in order, i.e., {b1 b2 b3})
% OUTPUTS:
%   combined [FT struct] - concatenated data structs

%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Concatenate Blocks
% Check properties of blocks array
same_chan = zeros([1 numel(blocks)-1]);
same_fs   = zeros([1 numel(blocks)-1]);
one_trl   = zeros([1 numel(blocks)]);
for b_ix = 1:numel(blocks)-1
    % Check if channel labels match
    if isempty(setdiff(blocks{b_ix+1}.label,blocks{b_ix}.label)) && isempty(setdiff(blocks{b_ix}.label,blocks{b_ix+1}.label))
        same_chan(b_ix) = 1;
    end
    % Check if same sampling rate
    if blocks{b_ix}.fsample==blocks{b_ix+1}.fsample
        same_fs(b_ix) = 1;
    end
    % Check for nested blocks
    if numel(blocks{b_ix})==1
        one_trl(b_ix) = 1;
    end
end
if numel(blocks{end})==1
    one_trl(end) = 1;
end
assert(iscell(blocks), 'blocks input is not a cell');
assert(all(same_chan), 'blocks have different channels!');
assert(all(same_fs), 'blocks have different sample rates!');
assert(all(one_trl), 'at least one block has more than one trial!');

% If only one block, just return it
if numel(blocks)<2
    combined = blocks{1};
    fprintf('===================================================================================\n');
    fprintf('WARNING: Only one block submitted to fn_concat_blocks, returning with no change!!!\n');
    fprintf('===================================================================================\n');
    
% If 2 or more blocks, combine them    
else
    % Initialize matrices
    n_blocks   = numel(blocks);
    block_lens = zeros([1 n_blocks]);
    end_times  = zeros([1 n_blocks]);
    for b_ix = 1:n_blocks
        end_times(b_ix) = blocks{b_ix}.time{1}(end);
        block_lens(b_ix) = numel(blocks{b_ix}.time{1});
    end
    time_step = 1/blocks{1}.fsample;
    
    data_concat = zeros([numel(blocks{1}.label) sum(block_lens)]);
    time_concat = zeros([1 sum(block_lens)]);
    
    % Combine data and time vectors
    data_concat(:,1:block_lens(1)) = blocks{1}.trial{1};
    time_concat(1,1:block_lens(1)) = blocks{1}.time{1};
    for b_ix = 2:n_blocks
        data_concat(:,sum(block_lens(1:b_ix-1))+1:sum(block_lens(1:b_ix-1))+block_lens(b_ix)) = blocks{b_ix}.trial{1};
        time_concat(:,sum(block_lens(1:b_ix-1))+1:sum(block_lens(1:b_ix-1))+block_lens(b_ix)) = blocks{b_ix}.time{1}+sum(end_times(1:b_ix-1))+time_step;
    end
    
    % Create Fieldtrip struct
    combined = blocks{1};
    combined.trial = {data_concat};
    combined.time = {time_concat};
    combined.sampleinfo = [1 size(data_concat,2)];
end

end
