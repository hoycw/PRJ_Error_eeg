function [labels, diff_pairs, colors, line_styles] = fn_condition_diff_label_styles(factor_name)
%% Converts a set of conditions into difference pairs (labels, plotting colors/styles)

switch factor_name
    case 'Holroyd'
        labels   = {'HL-EW', 'EL-HW'};
        [cond_lab,~,~,~] = fn_condition_label_styles(factor_name);
        diff_pairs = {[find(strcmp(cond_lab,'HdLs')) find(strcmp(cond_lab,'EzWn'))]...
            [find(strcmp(cond_lab,'EzLs')) find(strcmp(cond_lab,'HdWn'))]};
        colors = {[31 120 180]./256, [227 26 28]./256};
        line_styles = {'-', '-'};
    case 'DifOutUE'
        labels   = {'HW-EW', 'EL-HL'};
        [cond_lab,~,~,~] = fn_condition_label_styles(factor_name);
        diff_pairs = {[find(strcmp(cond_lab,'HdWn')) find(strcmp(cond_lab,'EzWn'))]...
            [find(strcmp(cond_lab,'EzLs')) find(strcmp(cond_lab,'HdLs'))]};
        colors = {[31 120 180]./256, [227 26 28]./256};
        line_styles = {'-', '-'};
    case 'DifOutdO'
        labels   = {'EW-HW', 'EL-HL'};
        [cond_lab,~,~,~] = fn_condition_label_styles(factor_name);
        diff_pairs = {[find(strcmp(cond_lab,'EzWn')) find(strcmp(cond_lab,'HdWn'))]...
            [find(strcmp(cond_lab,'EzLs')) find(strcmp(cond_lab,'HdLs'))]};
        colors = {[31 120 180]./256, [227 26 28]./256};
        line_styles = {'-', '-'};
    case 'DifOutWL'
        labels   = {'EW-EL', 'HW-HL'};
        [cond_lab,~,~,~] = fn_condition_label_styles(factor_name);
        diff_pairs = {[find(strcmp(cond_lab,'EzWn')) find(strcmp(cond_lab,'EzLs'))]...
            [find(strcmp(cond_lab,'HdWn')) find(strcmp(cond_lab,'HdLs'))]};
        colors = {[31 120 180]./256, [227 26 28]./256};
        line_styles = {'-', '-'};
    case 'DifOutS'
        labels   = {'HW-EW', 'EL-HL', 'Su'};
        [cond_lab,~,~,~] = fn_condition_label_styles(factor_name);
        diff_pairs = {[find(strcmp(cond_lab,'HdWn')) find(strcmp(cond_lab,'EzWn'))]...
            [find(strcmp(cond_lab,'EzLs')) find(strcmp(cond_lab,'HdLs'))]...
            [find(strcmp(cond_lab,'Su')) 0]};
        colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        line_styles = {'-', '-', '-'};
    case 'DifOutSx'
        labels   = {'HW-EW', 'EL-HL', 'Su-EWHL'};
        [cond_lab,~,~,~] = fn_condition_label_styles(factor_name);
        diff_pairs = {[find(strcmp(cond_lab,'HdWn')) find(strcmp(cond_lab,'EzWn'))]...
            [find(strcmp(cond_lab,'EzLs')) find(strcmp(cond_lab,'HdLs'))]...
            {find(strcmp(cond_lab,'Su')) {[find(strcmp(cond_lab,'EzWn')) find(strcmp(cond_lab,'HdLs'))]}}};
        colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        line_styles = {'-', '-', '-'};
    case 'DiffOutStdTar'
        labels  = {'Tar-Std'}
        [cond_lab,~,~,~] = fn_condition_label_styles('Odd')
        diff_pairs= {[find(strcmp(cond_lab,'tar')) find(strcmp(cond_lab,'std'))]};
        colors = {[31 120 180]./256};
        line_styles = {'-'};
    otherwise
        error(['Difference wave not defined for factor ' factor_name]);
end

end