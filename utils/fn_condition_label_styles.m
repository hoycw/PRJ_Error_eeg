function [labels, colors, line_styles, markers] = fn_condition_label_styles(factor_name)
%% Converts the name of a set of conditions into labels, plotting colors/styles
% condition_name: [str] 'EH'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3
%   light red: [251 154 153]
%   dark red: [227 26 28]
%   light blue: [166 206 227]
%   dark blue: [31 120 180]
%   light green: [178 223 138]
%   dark green: [51 160 44]

% if length(cond_lab) == 1
switch factor_name
    case 'Odd'
        labels = {'std','tar','odd'};
        colors = {[168 180 165]./256, [144 205 229]./256, [142 82 126]./256};
        line_styles = {'-', '-', '-'};
        markers = {'o', 'o', 'o'};
    case 'DifOdd_OS'
        error('fix color label mismatch');
        labels = {'odd', 'std'};
        colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        line_styles = {'-', '-'};
        markers = {'o', 'o'};
    case 'DifOdd_TS'
        error('fix color label mismatch');
        labels = {'tar','std'};
        colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        line_styles = {'-', '-'};
        markers = {'o', 'o'};
    case 'DifOdd_OT'
        error('fix color label mismatch');
        labels = {'odd', 'tar'};
        colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        line_styles = {'-', '-'};
        markers = {'o', 'o'};
    case 'Dif'
        labels = {'Ez', 'Hd'};
        colors = {[228,26,28]./256, [55,126,184]./256};
        line_styles = {'-', '-'};    % colors for cond_lab plotting
        markers = {'o', 'd'};
    case 'Out'
        labels = {'Wn', 'Ls'};
        colors = {[31 120 180]./256, [227 26 28]./256};
        line_styles = {'-', '-'};    % colors for cond_lab plotting
        markers = {'o', 'o'};
    case {'OutS','FB'}
        labels = {'Wn', 'Ls', 'Su'};
        colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
        markers = {'o', 'o', 'o'};
    case {'DifOut','DifOutUE','DifOutWL','DifOutdO','DifOutDO','Holroyd'}
        labels = {'EzWn', 'EzLs', 'HdWn', 'HdLs'};
        colors = {[166 206 227]./256, [251 154 153]./256, [31 120 180]./256, [227 26 28]./256};
        line_styles = {'-', '-', '--', '--'};
        markers = {'o', 'o', 'd', 'd'};
    case {'DifOutS', 'DifOutSx'}
        labels = {'EzWn', 'EzLs', 'HdWn', 'HdLs', 'Su'};
        colors = {[166 206 227]./256, [251 154 153]./256, ...
                  [31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        line_styles = {'-', '-', '--', '--', ':'};
        markers = {'o', 'o', 'd', 'd', '*'};
    case {'DifOutSur','DifFB'}
        labels = {'EzWn', 'EzLs', 'EzSu', 'HdWn', 'HdLs', 'HdSu'};
        colors = {[166 206 227]./256, [251 154 153]./256, [178 223 138]./256, ...
                  [31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        line_styles = {'-', '-','-', '--', '--','--'};
        markers = {'o', 'o', 'o', 'd', 'd', 'd'};
    case 'Tim'
        labels = {'Er', 'Lt'};
        colors = {[247,104,161]./256, [122,1,119]./256};    % pink and purple
        line_styles = {'-', '-'};
        markers = {'o', 'o'};
    case 'EzOut'
        labels = {'EzWn', 'EzLs'};
        colors = {[166 206 227]./256, [251 154 153]./256};
        line_styles = {'-', '-'};
        markers = {'o', 'o'};
    case 'EzOutS'
        labels = {'EzWn', 'EzLs', 'EzSu'};
        colors = {[166 206 227]./256, [251 154 153]./256, [178 223 138]./256};
        line_styles = {'-', '-', '-'};
        markers = {'o', 'o', 'o'};
    case 'HdOut'
        labels = {'HdWn', 'HdLs'};
        colors = {[31 120 180]./256, [227 26 28]./256};
        line_styles = {'-', '-'};
        markers = {'d', 'd'};
    case 'HdOutS'
        labels = {'HdWn', 'HdLs', 'HdSu'};
        colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        line_styles = {'-', '-', '-'};
        markers = {'d', 'd', 'd'};
    otherwise
        error(strcat('Only one, unrecognized condition offered: ',factor_name));
end

end
