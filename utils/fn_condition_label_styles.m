function [labels, colors, line_styles, markers] = fn_condition_label_styles(factor_name)
%% Converts the name of a set of conditions into labels, plotting colors/styles
% condition_name: [str] 'EH'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

% if length(cond_lab) == 1
switch factor_name
%    case 'CNI'
%        labels = {'con', 'neu', 'inc'};
%        colors = {[55,126,184]./256, [0 0 0], [228,26,28]./256};
%        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
    case 'Dif'
        labels = {'Ez', 'Hd'};
        colors = {[228,26,28]./256, [55,126,184]./256};
        line_styles = {'-', '-'};    % colors for cond_lab plotting
        markers = {'o', 'd'};
    case 'Out'
        labels = {'Wn', 'Ls'};
        colors = {[55,126,184]./256, [228,26,28]./256};
        line_styles = {'-', '-'};    % colors for cond_lab plotting
        markers = {'o', 'o'};
    case 'DifOut'
        labels = {'EzWn', 'EzLs', 'HdWn', 'HdLs'};
        colors = {[55,126,184]./256, [228,26,28]./256, [55,126,184]./256, [228,26,28]./256};
        line_styles = {'-', '-', '-', '-'};
        markers = {'o', 'o', 'd', 'd'};
    case 'Tim'
        labels = {'Er', 'Lt'};
        colors = {[1 1 1], [0 0 0]};
        line_styles = {'-', '-'};
        markers = {'o', 'o'};
    case 'EzOut'
        labels = {'EzWn', 'EzLs'};
        colors = {[55,126,184]./256, [228,26,28]./256};
        line_styles = {'-', '-'};
        markers = {'o', 'o'};
    case 'HdOut'
        labels = {'HdWn', 'HdLs'};
        colors = {[55,126,184]./256, [228,26,28]./256};
        line_styles = {'-', '-'};
        markers = {'d', 'd'};
%    case 'pcon'
%        labels = {'mcon', 'same', 'minc'};
%        colors = {[55,126,184]./256, [0 0 0], [228,26,28]./256};
%        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
%    case 'pcon_CI'
%        labels = {'con_mcon', 'con_minc', 'inc_mcon', 'inc_minc'};
%        colors = {[0 0 0], [0 0 0], [228,26,28]./256, [228,26,28]./256};    % colors for cond_lab plotting
%        line_styles = {'-', '--', '-','--'};    % colors for cond_lab plotting
%     case 'conseq'
%         cond_id = 'conseq';
    otherwise
        error(strcat('Only one, unrecognized condition offered: ',factor_name));
end
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

end
