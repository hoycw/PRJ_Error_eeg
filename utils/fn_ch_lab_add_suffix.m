function [suffix_labels] = fn_ch_lab_add_suffix(labels,suffix)

if isempty(labels)
    suffix_labels = labels;
end

% Add '-' to each label
for i = 1:numel(labels);
    suffix_labels{i} = [labels{i} suffix];
end

end
