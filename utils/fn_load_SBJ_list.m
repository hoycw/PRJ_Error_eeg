function SBJs = fn_load_SBJ_list(SBJ_id)
%% Load SBJ list from file
% INPUTS:
%   SBJ_id [str] - label of the SBJ list
% OUTPUTS:
%   SBJs [cell array] - string list of SBJs

[root_dir, ~] = fn_get_root_dir();

% Select SBJs
sbj_file = fopen([root_dir 'PRJ_Error_EEG/scripts/SBJ_lists/' SBJ_id '.sbj']);
tmp = textscan(sbj_file,'%s');
fclose(sbj_file);
SBJs = tmp{1};

end
