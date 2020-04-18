function [labels, names, colors, line_styles] = fn_regressor_label_styles(model_id)
%% Converts the name of a model into regressor labels, plotting colors/styles
% No performance metrics, but no includes categorical outcomes
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

%% List of possible regressors and their colors
regressors  = {...
    'RVal','RMag','RLik',... % Categorical reward features
    'pWin','bAcc','rAcc','rAcc10','score',... % Outcome predictors
    'sRPE','uRPE',... % Reward Feedback
    'OLik',... % Outcome Likelihood
    'ITI', ... % preceding inter-trial interval
    'SBJ' ... % SBJ only null model
    };
regressor_names = {...
    'Reward Valence','Reward Magnitude','Reward Likelihood',...
    'Prob(Win)','Block Accuracy','Rolling Accuracy','Rolling Accuracy 10','Total Score',...
    '+/- Pred Err','abs Pred Err',...
    'Outcome Likelihood',...
    'ITI', ...
    'SBJ' ...
    };
% Colors with just RL model:
regressor_colors = {...
    [209 151 105]./255, [118 160 156]./255, [152 78 163]./255,... % tan, teal, purple
    [0 0 0], [0 0 0], [0 0 0], [0 0 0], [0.4 0.4 0.4],... % dark blue, blacks, gray
    [209 151 105]./255, [118 160 156]./255, ... % tan, teal
    [152 78 163]./255,... % purple
    [0.3 0.3 0.3], ...   % dark gray
    [0.3 0.3 0.3] ...   % dark gray
    };

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
    % Categorical Outcome Features
    case 'RV'
        labels = {'RVal'};
    case 'RVL'
        labels = {'RVal','RLik'};
    case 'RVM'
        labels = {'RVal','RMag'};
    case 'RVLM'
        labels = {'RVal','RLik','RMag'};
    
    % RL models with pWin
    case 'RPE'
        labels = {'pWin','sRPE','uRPE'};
    case 'sRPE'
        labels = {'pWin','sRPE'};
    case 'uRPE'
        labels = {'pWin','uRPE'};
    case 'RPEOL'
        labels = {'pWin','sRPE','uRPE','OLik'};
    case 'RPEiti'
        labels = {'pWin','sRPE','uRPE','ITI'};
    case 'RPEscr'
        labels = {'pWin','score','sRPE','uRPE'};
    
    % RL models with alternative accuracy
%     case 'RLbA'
%         labels = {'bAcc','sRPE','uRPE'};
%     case 'RLbApW'
%         labels = {'bAcc','pWin','sRPE','uRPE'};
%     case 'RLrA10'
%         labels = {'rAcc10','sRPE','uRPE'};
    case 'RLrA'
        labels = {'rAcc','sRPE','uRPE'};
%     case 'RLrApW'
%         labels = {'rAcc','pWin','sRPE','uRPE'};
    
    % Null SBJ only control model
    case 'SBJonly'
        labels = {'SBJ'};
        
    % Previous Trial Regressors:
    case 'RLpT'
        error('why running with previous trial regressors?');
        labels = {'pWin','sRPE','uRPE','sTaEr','psTaEr'};
    case 'RLpRTlD'
        error('why running with previous trial regressors?');
        labels = {'pWin','sRPE','uRPE','sTaEr','psTaEr','sThPr'};
    case 'RLpRTulD'
        error('why running with previous trial regressors?');
        labels = {'pWin','sRPE','uRPE','sTaEr','psTaEr','sThPr','uThPr'};
%     case 'RLpRTiD'
%         error('why running with previous trial regressors?');
%         labels = {'pWin','sRPE','uRPE','sTaEr','psTaEr','iThr'};
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
    if any(strcmp(labels{reg_ix},'psTaEr'))%{'psTaEr','uThPr'} % strcmp(labels{reg_ix}(1),'u') || 
        line_styles{reg_ix} = ':';
    elseif strcmp(labels{reg_ix},'p2sTaEr')
        line_styles{reg_ix} = '-.';
    else
        line_styles{reg_ix} = '-';
    end
end

end