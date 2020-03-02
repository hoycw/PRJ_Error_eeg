function [labels, names, colors, line_styles] = fn_regressor_label_styles(model_id)
%% Converts the name of a model into regressor labels, plotting colors/styles
% Future:
%   add unsigned distance metrics
%   add categorical 0/1 if becomes necessary
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

%% List of possible regressors and their colors
regressors  = {...
    'pWin','sPE','uPE',... % Reward Feedback
    'sTarD','psTarD','p2sTarD',... % Performance (signed)
    'uTarD','puTarD','p2uTarD',... % Performance (unsigned)
    'sThrD','uThrD' ... %'iThrD', ... % Performance (threshold)
    };
regressor_names = {...
    'Prob(Win)','+/- Pred Err','abs Pred Err',...
    '+/- Target Dist','+/- Target Dist (n-1)','+/- Target Dist (n-2)',...
    'abs Target Dist','abs Target Dist (n-1)','abs Target Dist (n-2)',...
    '+/- Thresh Dist','abs Thresh Dist', ...%'+/- Thresh Dist (^-1)'...
    };
regressor_colors = {...
    [0 0 0], [118 160 156]./255, [209 151 105]./255, ... % black, teal, tan
    [152 78 163]./255, [152 78 163]./255, [152 78 163]./255, ...% purple
    [247 129 191]./255, [247 129 191]./255, [247 129 191]./255, ...% pink
    [97 61 46]./255, [189 65 45]./255 ...%brown, red
    };

%   Original Plots: {pWin, sPE, uPE, tRT (3x), lDist (3x)}
% regressor_colors = {...
%     [152 78 163]./255, [255 127 0]./255, [166 86 40]./255,...%purple, orange, brown
%     [247 129 191]./255, [247 129 191]./255, [247 129 191]./255, ...%pink, pink, pink
%     [0.5 0.5 0.5], [0.5 0.5 0.5], [0.5 0.5 0.5] ...%gray, gray, gray
%     };

% Color options:
% Taken from: colorbrewer2.org, qualitative, 5-class Set1
%   purple: [152 78 163]
%   orange: [255 127 0]
%   brown:  [166 86 40]
%   yellow: [255 255 51]
%   pink:   [247 129 191]
%   red:    [228 26 28]
%   green:  [55 126 184]
%   blue:   [77 175 74]
%   gray:   [0.5 0.5 0.5]
% Found in Online Image "color_combinations" on Desktop
%   Antiquated Aqua:    [118 160 156]
%   Warm Sunglow (tan): [209 151 105]
%   Tarrazzo Brown:     [97 61 46]
%   Tomato Tango:       [189 65 45]

%% Convert model_id into set of conditions
switch model_id
    case 'RL'
        labels = {'pWin','sPE','uPE'};
    case 'RLsTarD'
        labels = {'pWin','sPE','uPE','sTarD'};
    case 'RLTarD'
        labels = {'pWin','sPE','uPE','sTarD','uTarD'};
    case 'RLThrD'
        labels = {'pWin','sPE','uPE','sThrD','uThrD'};
    case 'RLfullD'
        labels = {'pWin','sPE','uPE','sTarD','uTarD','sThrD','uThrD'};
    case 'RL3D'
        labels = {'pWin','sPE','uPE','sTarD','uTarD','uThrD'};
    case 'RLsD'
        labels = {'pWin','sPE','uPE','sTarD','sThrD'};
    case 'RLuD'
        labels = {'pWin','sPE','uPE','uTarD','uThrD'};
        
    % Pre-Feedback Models
    case 'pWTarD'
        labels = {'pWin','sTarD','uTarD'};
    case 'pWallD'
        labels = {'pWin','sTarD','uTarD','sThrD','uThrD'};
        
    % Previous Trial Regressors:
    case 'RLpT'
        error('why running with previous trial regressors?');
        labels = {'pWin','sPE','uPE','sTarD','psTarD'};
    case 'RLpRTlD'
        error('why running with previous trial regressors?');
        labels = {'pWin','sPE','uPE','sTarD','psTarD','sThrD'};
    case 'RLpRTulD'
        error('why running with previous trial regressors?');
        labels = {'pWin','sPE','uPE','sTarD','psTarD','sThrD','uThrD'};
%     case 'RLpRTiD'
%         error('why running with previous trial regressors?');
%         labels = {'pWin','sPE','uPE','sTarD','psTarD','iThrD'};
    otherwise
        error(strcat('Unknown model_id: ',model_id));
end

%% Assign colors and line styles
names       = cell(size(labels));
colors      = cell(size(labels));
line_styles = cell(size(labels));
for reg_ix = 1:numel(labels)
    names{reg_ix}  = regressor_names{strcmp(labels{reg_ix},regressors)};
    colors{reg_ix} = regressor_colors{strcmp(labels{reg_ix},regressors)};
    if any(strcmp(labels{reg_ix},'psTarD'))%{'psTarD','uThrD'}
        line_styles{reg_ix} = ':';
    elseif strcmp(labels{reg_ix},'p2sTarD')
        line_styles{reg_ix} = '-.';
    else
        line_styles{reg_ix} = '-';
    end
end

end
