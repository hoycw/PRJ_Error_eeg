function [prefix_labels] = fn_ch_lab_add_prefix(labels,prefix)

if isempty(labels)
    prefix_labels = labels;
end

% Add '-' to each label
for i = 1:numel(labels);
    prefix_labels{i} = [prefix labels{i}];
end

end
