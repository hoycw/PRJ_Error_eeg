function [grp_labels, colors, line_styles] = fn_group_label_styles(model_id)
%% Converts the name of a set of conditions into labels, plotting colors/styles
% condition_name: [str] 'CNI', 'CI', 'pcon'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

%% List of possible labels and their colors
factors  = {'Dif','Out','Tim','Dif*Out','RT'};
fact_colors = {[55 126 184]./256, [228 26 28]./256, [77 175 74]./256, [152 78  163]./256, [0.5 0.5 0.5]};
% Taken from: colorbrewer2.org, qualitative, 5-class Set1
%   currently: green, red, blue, purple, gray
%   orange for later: [255 127 0]

%% Convert model_id into set of conditions
switch model_id
    case 'DifOutTim'
        grp_labels = {'Dif','Out','Tim'};
%     case 'rRT_CNI_pcon'
%         labels = {'CNI', 'pcon'};
    case 'corrRT_DifOutTim'
        grp_labels = {'Dif','Out','Tim'};
    case 'corrRT_DifOutTimDO'
        grp_labels = {'Dif','Out','Tim','Dif*Out'};
    case 'DifOut'
        grp_labels = {'Dif','Out'};
    case 'corrRT_DifOut'
        grp_labels = {'Dif','Out'};
    case 'RT'
        grp_labels = {'RT'};
    otherwise
        error(strcat('Unknown model_id: ',model_id));
end

% Assign colors and line styles
colors = {};
line_styles = {};
for cond_ix = 1:numel(grp_labels)
    colors{cond_ix} = fact_colors{strmatch(grp_labels{cond_ix},factors,'exact')};
    if isempty(strfind(grp_labels{cond_ix},'*'))
        line_styles{cond_ix} = '-';     % Main effects
    else
        line_styles{cond_ix} = '--';    % Interactions
    end
end

end

%% Extra conditions
%     case 'CI'
%         labels = {'con', 'inc'};
%         colors = {[55,126,184]./256, [228,26,28]./256};
%         line_styles = {'-', '-'};    % colors for cond_lab plotting
%     case 'pcon'
%         labels = {'mcon', 'same', 'minc'};
%         colors = {[55,126,184]./256, [0 0 0], [228,26,28]./256};
%         line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
%     case 'pcon_CI'
%         labels = {'con_mcon', 'con_minc', 'inc_mcon', 'inc_minc'};
%         colors = {[0 0 0], [0 0 0], [228,26,28]./256, [228,26,28]./256};    % colors for cond_lab plotting
%         line_styles = {'-', '--', '-','--'};    % colors for cond_lab plotting
%     case 'conseq'
%         cond_id = 'conseq';

% else
%     cond_id = 'cst';
%     for c_ix = 1:length(cond_lab)
%         cond_id = [cond_id fn_convert_condition_lab2num(cond_lab{c_ix})];
%     end
%     cond_colors = manual_cond_colors;
%     if length(cond_lab)~=length(cond_colors)
%         error('Mismatched condition labels and colors');
%     end
% end
