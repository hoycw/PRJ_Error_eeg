function [grp_labels, colors, line_styles] = fn_group_label_styles(model_id)
%% Converts the name of a set of conditions into labels, plotting colors/styles
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

%% List of possible labels and their colors
factors  = {'Dif','Out','Tim','Dif*Out','Sur'};
fact_colors = {[152 78 163]./255, [255 127 0]./255, [255 255 51]./255, [166 86 40]./255, [166 86 40]./255};
% Newer (different than RGB for 3 feedback conditions:
%   purple, orange, yellow, brown (brown again for Sur/DO repeat)
%   pink if needed: [247 129 191]
% Taken from: colorbrewer2.org, qualitative, 5-class Set1
% OLD:
%   fact_colors = {[55 126 184]./256, [228 26 28]./256, [77 175 74]./256, [152 78  163]./256, [0.5 0.5 0.5]};
%       green, red, blue, purple, gray
%   orange for later: [255 127 0]

%% Convert model_id into set of conditions
switch model_id
    case 'DifOutTim'
        grp_labels = {'Dif','Out','Tim'};
    case 'corrRT_DifOutTim'
        grp_labels = {'Dif','Out','Tim'};
    case 'corrRT_DifOutTimDO'
        grp_labels = {'Dif','Out','Tim','Dif*Out'};
    case 'DifOut'
        grp_labels = {'Dif','Out'};
    case 'DifOutSur'
        grp_labels = {'Dif','Out','Sur'};
    case 'DifOutS'
        grp_labels = {'Dif','Out'};
    case 'corrRT_DifOut'
        grp_labels = {'Dif','Out'};
    case 'Out'
        grp_labels = {'Out'};
    case 'OutS'
        grp_labels = {'OutS'};
    case 'Dif'
        grp_labels = {'Dif'};
    case 'RT'
        grp_labels = {'RT'};
    otherwise
        error(strcat('Unknown model_id: ',model_id));
end

% Assign colors and line styles
colors = cell(size(grp_labels));
line_styles = cell(size(grp_labels));
for cond_ix = 1:numel(grp_labels)
    colors{cond_ix} = fact_colors{strcmp(grp_labels{cond_ix},factors)};
    if isempty(strfind(grp_labels{cond_ix},'*'))
        line_styles{cond_ix} = '-';     % Main effects
    else
        line_styles{cond_ix} = '--';    % Interactions
    end
end

end
