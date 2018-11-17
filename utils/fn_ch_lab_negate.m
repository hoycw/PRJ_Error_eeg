function [negative_labels] = fn_ch_lab_negate(labels)

if isempty(labels)
    negative_labels = labels;
end

% Add '-' to each label
for i = 1:numel(labels);
    negative_labels{i} = ['-' labels{i}];
end

end
