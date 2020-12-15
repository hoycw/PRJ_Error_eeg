function [grp_labels, grp_names, colors, line_styles] = fn_group_label_styles(model_id)
%% Converts the name of a model into a set of group labels, plotting colors/styles
%   These are supersets of specific conditions, akin to factors in an ANOVA
%       difficulty (Dif): easy, hard
%       outcome (Out): win, neutral, loss
%       timing (Tim): early, late
%       tolerance (Tol)
%       reaction time (RT)
%       interactions indicated by A*B (e.g., Dif*Out)
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

%% List of possible labels and their colors
% Factor codes
factors  = {'OB','Dif','Tol','Out','OutS','FB','Tim','Dif*Out','Tol*Out','Sur'};

% Factor labels (longer strings for clear plotting legends)
factor_names = {'Oddball','Difficulty', 'Tolerance', 'Outcome', 'Outcome Surprise', 'Feedback', ...
    'Timing', 'Difficulty*Outcome', 'Tolerance*Outcome', 'Surprise'};

% Factor colors ([R G B] scaled to 0-1)
factor_colors = {[0.3 0.3 0.3], [152 78 163]./255, [152 78 163]./255, [255 127 0]./255, [255 127 0]./255, [255 127 0]./255,...
    [255 255 51]./255, [166 86 40]./255, [166 86 40]./255, [166 86 40]./255};
% Newer (different than RGB for 3 feedback conditions:
%   purple, orange, yellow, brown (brown again for Sur/DO repeat)
%   pink if needed: [247 129 191]
% Taken from: colorbrewer2.org, qualitative, 5-class Set1

%% Convert model_id into set of conditions
switch model_id
    case 'OB'
        grp_labels = {'OB'};
        
    case 'Dif'
        grp_labels = {'Dif'};
    case 'Out'
        grp_labels = {'Out'};
    case 'DifFB'
        grp_labels = {'Dif','FB'};
    case {'DifOut', 'corrRT_DifOut'}
        grp_labels = {'Dif','Out'};
        
    % OLDER unused factors:
    case {'DifOutTim', 'corrRT_DifOutTim'}
        grp_labels = {'Dif','Out','Tim'};
    case 'corrRT_DifOutTimDO'
        grp_labels = {'Dif','Out','Tim','Dif*Out'};
    case 'DifOutDO'
        grp_labels = {'Dif','Out','Dif*Out'};
    case 'TolOut'
        grp_labels = {'Tol','Out'};
    case 'TolOutDO'
        grp_labels = {'Tol','Out','Tol*Out'};
    case 'DifOutSur'
        grp_labels = {'Dif','Out','Sur'};
    case 'RT'
        grp_labels = {'RT'};
    case 'DiffOddOS'
        grp_labels = {'OddOS'};
    % 'OutS' should now be called 'FB'
%     case 'DifOutS'
%         grp_labels = {'Dif','OutS'};
%     case 'OutS'
%         grp_labels = {'OutS'};
    
    otherwise
        error(strcat('Unknown model_id: ',model_id));
end

%% Assign colors and line styles
grp_names   = cell(size(grp_labels));
colors      = cell(size(grp_labels));
line_styles = cell(size(grp_labels));
for grp_ix = 1:numel(grp_labels)
    grp_names{grp_ix} = factor_names{strcmp(grp_labels{grp_ix},factors)};
    colors{grp_ix}    = factor_colors{strcmp(grp_labels{grp_ix},factors)};
    % Define Line Styles
    if isempty(strfind(grp_labels{grp_ix},'*'))
        line_styles{grp_ix} = '-';     % Main effects
    else
        line_styles{grp_ix} = '--';    % Interactions
    end
end

end
