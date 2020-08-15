function [suffix_labels] = fn_ch_lab_add_suffix(labels,suffix)
% Add suffix string to end of channel labels
% INPUTS:
%   labels [cell array] - string labels of channels
%   suffix [str] - string to append to end of each channel label
% OUTPUTS:
%   suffix_labels [cell array] - channel labels with suffix appended

if isempty(labels)
    suffix_labels = labels;
end

% Add suffix to each label
for i = 1:numel(labels);
    suffix_labels{i} = [labels{i} suffix];
end

end
