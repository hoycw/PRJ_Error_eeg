function [labels, names, colors, line_styles, markers] = fn_condition_label_styles(grp_id)
%% Converts the name of a group of conditions into labels, plotting colors/styles
% condition_name: [str] 'EH'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3
%   light red: [251 154 153]
%   dark red: [227 26 28]
%   light blue: [166 206 227]
%   dark blue: [31 120 180]
%   light green: [178 223 138]
%   dark green: [51 160 44]

%% List of possible labels and their colors
conditions  = {...
    'Std','Tar','Odd',...
    'Ez','Hd','Wn','Ls','Su',...
    'EzWn','EzLs','EzSu','HdWn','HdLs','HdSu',...
    'All','Pos','Neg',...
    'Er','Lt','ErQ1','ErQ2','MdQ3','LtQ4','LtQ5'};
condition_names = {...
    'Standard','Target','Oddball',...
    'Easy','Hard','Win','Loss','Surprise',...
    'Easy Win','Easy Loss','Easy Surprise','Hard Win','Hard Loss','Hard Surprise',...
    'All Trials', 'Positive', 'Negative',...
    'Early','Late','Early Q1','Early Q2','Middle Q3','Late Q4','Late Q5'};
condition_colors = {...
    [168 180 165]./256, [144 205 229]./256, [142 82 126]./256, ...
    [228,26,28]./256, [55,126,184]./256, [31 120 180]./256, [227 26 28]./256, [51 160 44]./256, ...
    [166 206 227]./256, [251 154 153]./256, [178 223 138]./256, ...
    [31 120 180]./256, [227 26 28]./256, [51 160 44]./256, ...
    [0 0 0], [31 120 180]./256, [227 26 28]./256,...
    [166 97 26]./256, [1 133 113]./256, ... % brown, aqua
    [166 97 26]./256, [223 194 125]./256, [0.4 0.4 0.4], [128 205 193]./256, [1 133 113]./256 ...    % [brown, tan, gray, teal, aqua]
    };

%% Convert grp_id into set of conditions
switch grp_id
    % Oddball Conditions
    case 'Odd'
        labels = {'Std','Tar','Odd'};
        %colors = {[168 180 165]./256, [144 205 229]./256, [142 82 126]./256};
        %line_styles = {'-', '-', '-'};
        %markers = {'o', 'o', 'o'};
    case 'DifOdd_OS'
        %error('fix color label mismatch');
        labels = {'Odd', 'Std'};
        %colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        %line_styles = {'-', '-'};
        %markers = {'o', 'o'};
    case 'DifOdd_TS'
        %error('fix color label mismatch');
        labels = {'Tar','Std'};
        %colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        %line_styles = {'-', '-'};
        %markers = {'o', 'o'};
    case 'DifOdd_OT'
        %error('fix color label mismatch');
        labels = {'Odd', 'Tar'};
        %colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        %line_styles = {'-', '-'};
        %markers = {'o', 'o'};
    
    % ---------------------------------------------------------------------
    % Target Time Conditions
    case 'All'
        labels = {'All'};
    case 'Dif'
        labels = {'Ez', 'Hd'};
        %colors = {[228,26,28]./256, [55,126,184]./256};
        %line_styles = {'-', '-'};    % colors for cond_lab plotting
        %markers = {'o', 'd'};
    case 'Out'
        labels = {'Wn', 'Ls'};
        %colors = {[31 120 180]./256, [227 26 28]./256};
        %line_styles = {'-', '-'};    % colors for cond_lab plotting
        %markers = {'o', 'o'};
    case 'Val'
        labels = {'Pos', 'Neg'};
    case {'OutS','FB'}
        labels = {'Wn', 'Ls', 'Su'};
        %colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        %line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
        %markers = {'o', 'o', 'o'};
    case {'DifOut','DifOutUE','DifOutWL','DifOutdO','DifOutDO','Holroyd'}
        labels = {'EzWn', 'EzLs', 'HdWn', 'HdLs'};
        %colors = {[166 206 227]./256, [251 154 153]./256, [31 120 180]./256, [227 26 28]./256};
        %line_styles = {'-', '-', '--', '--'};
        %markers = {'o', 'o', 'd', 'd'};
    case {'DifOutS', 'DifOutSx'}
        labels = {'EzWn', 'EzLs', 'HdWn', 'HdLs', 'Su'};
        %colors = {[166 206 227]./256, [251 154 153]./256, ...
        %          [31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        %line_styles = {'-', '-', '--', '--', ':'};
        %markers = {'o', 'o', 'd', 'd', '*'};
    case {'DifOutSur','DifFB'}
        labels = {'EzWn', 'EzLs', 'EzSu', 'HdWn', 'HdLs', 'HdSu'};
        %colors = {[166 206 227]./256, [251 154 153]./256, [178 223 138]./256, ...
        %          [31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        %line_styles = {'-', '-','-', '--', '--','--'};
        %markers = {'o', 'o', 'o', 'd', 'd', 'd'};
    case 'Tar2'
        labels = {'Er', 'Lt'};
%     case 'Tar4'
%         labels = {'ErQ1','ErQ2','LtQ3','LtQ4'};
    case 'Tar5'
        labels = {'ErQ1','ErQ2','MdQ3','LtQ4','LtQ5'};
        %colors = {[247,104,161]./256, [122,1,119]./256};    % pink and purple
        %line_styles = {'-', '-'};
        %markers = {'o', 'o'};
    case 'EzOut'
        labels = {'EzWn', 'EzLs'};
        %colors = {[166 206 227]./256, [251 154 153]./256};
        %line_styles = {'-', '-'};
        %markers = {'o', 'o'};
    case 'EzOutS'
        labels = {'EzWn', 'EzLs', 'EzSu'};
        %colors = {[166 206 227]./256, [251 154 153]./256, [178 223 138]./256};
        %line_styles = {'-', '-', '-'};
        %markers = {'o', 'o', 'o'};
    case 'HdOut'
        labels = {'HdWn', 'HdLs'};
        %colors = {[31 120 180]./256, [227 26 28]./256};
        %line_styles = {'-', '-'};
        %markers = {'d', 'd'};
    case 'HdOutS'
        labels = {'HdWn', 'HdLs', 'HdSu'};
        %colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
        %line_styles = {'-', '-', '-'};
        %markers = {'d', 'd', 'd'};
    otherwise
        error(strcat('Only one, unrecognized condition offered: ',grp_id));
end

%% Assign colors and line styles
names       = cell(size(labels));
colors      = cell(size(labels));
line_styles = cell(size(labels));
for cond_ix = 1:numel(labels)
    names{cond_ix}  = condition_names{strcmp(labels{cond_ix},conditions)};
    colors{cond_ix} = condition_colors{strcmp(labels{cond_ix},conditions)};
    % Define Line Styles
    if isempty(strfind(labels{cond_ix},'Hd'))
        line_styles{cond_ix} = '-';     % Easy
    else
        line_styles{cond_ix} = '--';    % Hard
    end
    % Define Marker Styles
    if strcmp(labels{cond_ix},'Su')
        markers{cond_ix} = '*';     % Surprise (difficulty agnostic)
    elseif isempty(strfind(labels{cond_ix},'Hd'))
        markers{cond_ix} = 'o';     % Easy
    else
        markers{cond_ix} = 'd';    % Hard
    end
end

end
