function [negative_labels] = fn_ch_lab_negate(labels)
% Add '-' string to beginning of channel labels to remove them using ft_selectdata
% INPUTS:
%   labels [cell array] - string labels of channels
% OUTPUTS:
%   negative_labels [cell array] - channel labels with '-' appended to front

if isempty(labels)
    negative_labels = labels;
end

% Add '-' to each label
for i = 1:numel(labels);
    negative_labels{i} = ['-' labels{i}];
end

end
