function [bhv] = fn_load_behav_csv(csv_fname)
%% Load trial info csv exported from python
%   Also converts from python 0-based to MATLAB 1-based indexing
% INPUTS:
%   csv_fname [str] - full file path to behavioral csv
% OUTPUTS:
%   bhv [struct] - contains many fields of info on the trials

% CSV format:
%   Total_Trial [int] - trial # from 1:n_trials*n_blocks (+training+examples)
%   Block [int] - block #
%   Condition [str] - 'easy' or 'hard'
%   Hit [int] - 0/1 if this trial was correct
%   RT [float] - reaction time on this trial
%   Timestamp [float] - time of outcome logging (since start of experiment)
%   Tolerance [float] - tolerance around target interval for this trial
%   Trial [int] - trial # within a block (e.g., 1:n/trials/block_n)
%   Score [int] - score obtained on this trial (+/- 100)
%   ITI [float] - time since timestamp n-1, will be 0 for first trial in block, training data
%   ITI Type [str] - ['short','medium', or 'long']
%       only have medium in cases with 3 (not 4) ITIs (e.g., IR57)


% Open file
fprintf('\tReading behavioral csv file: %s\n',csv_fname);
bhv_file = fopen(csv_fname);

% Read field names
py_fields  = {'Total_Trial', 'Block', 'Feedback', 'RT', 'Sound', 'Timestamp', 'Tolerance', 'Trial', 'Hit',  ...
              'Score', 'bad_fb', 'Condition', 'ITI', 'ITI type'};
new_fields = {'trl_n', 'blk', 'fb', 'rt', 'sound', 'time', 'tol', 'blk_trl_n', 'hit', ...
              'score', 'bad_fb', 'cond', 'ITI', 'ITI_type'};
bhv_fields = textscan(bhv_file, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s', 1, 'Delimiter', ',');

% OLD filed formats:
%py_fields  = {'Total_Trial', 'Block', 'Condition', 'Hit', 'RT', 'Timestamp', 'Tolerance', 'Trial', 'Score',  ...
              % 'ITI', 'ITI type'};
%new_fields = {'trl_n', 'blk', 'blk_trl_n', 'fb', 'rt', 'tol', 'time', 'hit', ...
              %'score', 'bad_fb', 'cond', 'ITI', 'ITI_type'};

% Check that loaded fields match expected py_fields
if numel(py_fields)~=numel(bhv_fields)
    error('Mismatched fields in behav csv than expected!');
end
for f = 1:numel(py_fields)
   if ~strcmp(py_fields{f},bhv_fields{f}{1})
      error(['Mismatched field in behav csv and expected: ' py_fields{f} ' vs. ' bhv_fields{f}{1}]);
    end
end

% Read data
%   orig version: formatspec = '%d%d%s%d%f%f%f%d%d%f%s';
bhv_data = textscan(bhv_file,'%d %d %s %f %s %f %f %d %d %d %s %s %f %f',...
                'Delimiter',',','HeaderLines',1);
%bhv_data = textscan(bhv_file,'%d %d %d %s %f %f %f %d %d %s %s %f %f',...
               % 'Delimiter',',','HeaderLines',1);
fclose(bhv_file);

% Get list of good trials
%   Not used here because all trials in log should be in recording
%   This is leftover from iEEG when sometimes trigger (photodiode) is lost
%   so that trial must be ignored in the behavioral log
n_trials = size(bhv_data{1},1);
fprintf('\t\tFound %d trials in log file\n', n_trials);
% if ~isempty(ignore_trials)
%     fprintf('\t\tIgnoring %d trials\n', numel(ignore_trials));
%     good_trials = setdiff(1:n_trials,ignore_trials);
% else
%     good_trials = 1:n_trials;
% end

% Add data to struct, convert field names with spaces to underscore
for f_ix = 1:numel(bhv_fields)
%     bhv_fields{f_ix} = strrep(bhv_fields{f_ix}{1},' ','_');
    if numel(bhv_data{f_ix})==n_trials
        bhv.(new_fields{f_ix}) = bhv_data{f_ix};%(good_trials);
    else
        bhv.(new_fields{f_ix}) = bhv_data{f_ix};
    end
end

% Convert from python to MATLAB indexing
bhv.trl_n     = bhv.trl_n+1;
bhv.blk       = bhv.blk+1;
bhv.blk_trl_n = bhv.blk_trl_n+1;

% Convert bad_fb from string to binary 0/1
str_bad_fb = bhv.bad_fb;
bhv.bad_fb = zeros(size(bhv.bad_fb));
bhv.bad_fb(strcmp(str_bad_fb,'True')) = 1;

end
