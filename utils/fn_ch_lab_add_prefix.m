function [prefix_labels] = fn_ch_lab_add_prefix(labels,prefix)
% Add prefix string to beginning of channel labels
%   NOTE: Use fn_ch_lab_negate to add '-' to labels to remove them using ft_selectdata
% INPUTS:
%   labels [cell array] - string labels of channels
%   prefix [str] - string to append to front of each channel label
% OUTPUTS:
%   prefix_labels [cell array] - channel labels with prefix appended

if isempty(labels)
    prefix_labels = labels;
end

% Add prefix to each label
for i = 1:numel(labels);
    prefix_labels{i} = [prefix labels{i}];
end

end
