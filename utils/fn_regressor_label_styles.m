function [labels, names, colors, line_styles] = fn_regressor_label_styles(model_id)
%% Converts the name of a model into regressor labels, plotting colors/styles
% Future:
%   add unsigned distance metrics
%   add categorical 0/1 if becomes necessary
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

%% List of possible regressors and their colors
regressors  = {...
    'pWin','bAcc','rAcc','rAcc10','score',... % Outcome predictors
    'sPE','uPE',... % Reward Feedback
    'sTaEr','psTaEr','p2sTaEr',... % Performance (signed)
    'sTaPr',... % Target precision (signed)
    'uTaEr','puTaEr','p2uTaEr',... % Performance (unsigned)
    'uTaPr',... % Target precision (unsigned)
    'sThEr','uThEr', ... % Threshold Error
    'sThPr','uThPr', ... %'iThr', ... % Threshold Precision
    'ITI', ... % preceding inter-trial interval
    'SBJ' ... % SBJ only null model
    };
regressor_names = {...
    'Prob(Win)','Block Accuracy','Rolling Accuracy','Rolling Accuracy 10','Total Score',...
    '+/- Pred Err','abs Pred Err',...
    '+/- Target Err','+/- Target Err (n-1)','+/- Target Err (n-2)',...
    '+/- Target Prec',...
    'abs Target Err','abs Target Err (n-1)','abs Target Err (n-2)',...
    'abs Target Prec',...
    '+/- Thresh Err','abs Thresh Err', ...
    '+/- Thresh Prec','abs Thresh Prec', ...%'+/- Thresh Prec (^-1)'...
    'ITI', ...
    'SBJ' ...
    };
regressor_colors = {...
    [0 0 0], [0 0 0], [0 0 0], [0 0 0], [0.4 0.4 0.4],... % blacks, gray
    [118 160 156]./255, [37 52 148]./255, ... % teal, dark blue
    [152 78 163]./255, [152 78 163]./255, [152 78 163]./255, ...% purple
    [247 104 161]./255,... % medium pink
    [247 129 191]./255, [247 129 191]./255, [247 129 191]./255, ...% pink
    [122 1 119]./255,... % dark purple
    [161 218 180]./255, [44 127 184]./255, ... % teal, medium blue
    [209 151 105]./255, [189 65 45]./255, ...%tan, red
    [0.3 0.3 0.3], ...   % dark gray
    [0.3 0.3 0.3] ...   % dark gray
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
    % Main Model
    case 'RL4D'
        labels = {'pWin','sPE','uPE','sTaEr','uTaEr','sThPr','uThPr'};
        
    % RL models with pWin
    case 'RL'
        labels = {'pWin','sPE','uPE'};
    case 'RLiti'
        labels = {'pWin','sPE','uPE','ITI'};
    case 'RLscr'
        labels = {'pWin','score','sPE','uPE'};
    case 'RLTaEr'
        labels = {'pWin','sPE','uPE','sTaEr','uTaEr'};
    case 'RLThPr'
        labels = {'pWin','sPE','uPE','sThPr','uThPr'};
    case 'RL3D'
        labels = {'pWin','sPE','uPE','sTaEr','uTaEr','uThPr'};
    case 'RLsD'
        labels = {'pWin','sPE','uPE','sTaEr','sThPr'};
    case 'RLuD'
        labels = {'pWin','sPE','uPE','uTaEr','uThPr'};
    case 'RLsTauTh'
        labels = {'pWin','sPE','uPE','sTaEr','uThPr'};
    case 'RLuTasTh'
        labels = {'pWin','sPE','uPE','uTaEr','sThPr'};
    
    % RL with Single Performance
    case 'RLsTaEr'
        labels = {'pWin','sPE','uPE','sTaEr'};
    case 'RLuTaEr'
        labels = {'pWin','sPE','uPE','uTaEr'};
    case 'RLsThPr'
        labels = {'pWin','sPE','uPE','sThPr'};
    case 'RLuThPr'
        labels = {'pWin','sPE','uPE','uThPr'};
    
    % RL models with alternative accuracy
%     case 'RLbA'
%         labels = {'bAcc','sPE','uPE'};
%     case 'RLbApW'
%         labels = {'bAcc','pWin','sPE','uPE'};
%     case 'RLbA3D'
%         labels = {'bAcc','sPE','uPE','sTaEr','uTaEr','uThPr'};
%     case 'RLrA10'
%         labels = {'rAcc10','sPE','uPE'};
    case 'RLrA'
        labels = {'rAcc','sPE','uPE'};
%     case 'RLrApW'
%         labels = {'rAcc','pWin','sPE','uPE'};
%     case 'RLrA3D'
%         labels = {'rAcc','sPE','uPE','sTaEr','uTaEr','uThPr'};
    
    % Null SBJ only control model
    case 'SBJonly'
        labels = {'SBJ'};
        
    % Performance Models
    case 'pWTaEr'
        labels = {'pWin','sTaEr','uTaEr'};
    case 'pWThPr'
        labels = {'pWin','sThPr','uThPr'};
    case 'pW4D'
        labels = {'pWin','sTaEr','uTaEr','sThPr','uThPr'};
    case 'pWscr4D'
        labels = {'pWin','score','sTaEr','uTaEr','sThPr','uThPr'};
    case 'rATaEr'
        labels = {'rAcc','sTaEr','uTaEr'};
    case 'rAThPr'
        labels = {'rAcc','sThPr','uThPr'};
    case 'rA4D'
        labels = {'rAcc','sTaEr','uTaEr','sThPr','uThPr'};
    case 'rAscr4D'
        labels = {'rAcc','score','sTaEr','uTaEr','sThPr','uThPr'};
        
    % OLD JUNK:
%     case 'RLTa'
%         labels = {'pWin','sPE','uPE','sTaEr','uTaEr','sTaPr','uTaPr'};
%     case 'RLTaPr'
%         labels = {'pWin','sPE','uPE','sTaPr','uTaPr'};
%     case 'RLTh'
%         labels = {'pWin','sPE','uPE','sThEr','uThEr','sThPr','uThPr'};
%     case 'RLThEr'
%         labels = {'pWin','sPE','uPE','sThEr','uThEr'};
%     case 'RLallD'
%         labels = {'pWin','sPE','uPE','sTaEr','uTaEr','sTaPr','uTaPr','sThEr','uThEr','sThPr','uThPr'};
%     case 'RLmostD'
%         labels = {'pWin','sPE','uPE','sTaEr','uTaEr','sTaPr','uTaPr','uThEr','uThPr'};
    
    % Previous Trial Regressors:
    case 'RLpT'
        error('why running with previous trial regressors?');
        labels = {'pWin','sPE','uPE','sTaEr','psTaEr'};
    case 'RLpRTlD'
        error('why running with previous trial regressors?');
        labels = {'pWin','sPE','uPE','sTaEr','psTaEr','sThPr'};
    case 'RLpRTulD'
        error('why running with previous trial regressors?');
        labels = {'pWin','sPE','uPE','sTaEr','psTaEr','sThPr','uThPr'};
%     case 'RLpRTiD'
%         error('why running with previous trial regressors?');
%         labels = {'pWin','sPE','uPE','sTaEr','psTaEr','iThr'};
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

%% Regressors before adding error and precision:
% regressors  = {...
%     'pWin','bAcc','rAcc',... % Outcome predictors
%     'sPE','uPE',... % Reward Feedback
%     'sTaEr','psTaEr','p2sTaEr',... % Performance (signed)
%     'sTaPr',... % Target precision (signed)
%     'uTaEr','puTaEr','p2uTaEr',... % Performance (unsigned)
%     'uTaPr',... % Target precision (unsigned)
%     'sThPr','uThPr', ... %'iThr', ... % Performance (threshold)
%     'SBJ' ... % SBJ only null model
%     };
% regressor_names = {...
%     'Prob(Win)','Block Accuracy','Rolling Accuracy',...
%     '+/- Pred Err','abs Pred Err',...
%     '+/- Target Err','+/- Target Err (n-1)','+/- Target Err (n-2)',...
%     '+/- Target Prec',...
%     'abs Target Err','abs Target Err (n-1)','abs Target Err (n-2)',...
%     'abs Target Prec',...
%     '+/- Thresh Prec','abs Thresh Prec', ...%'+/- Thresh Prec (^-1)'...
%     'SBJ' ...
%     };
% regressor_colors = {...
%     [0 0 0], [0 0 0], [0 0 0], ... % black
%     [118 160 156]./255, [209 151 105]./255, ... % teal, tan
%     [152 78 163]./255, [152 78 163]./255, [152 78 163]./255, ...% purple
%     []./255,...
%     [247 129 191]./255, [247 129 191]./255, [247 129 191]./255, ...% pink
%     []./255,...
%     [97 61 46]./255, [189 65 45]./255, ...%brown, red
%     [0.2 0.2 0.2] ...   % dark gray
%     };
% both error/precision:
% regressor_colors = {...
%     [0 0 0], [0 0 0], [0 0 0], ... % black
%     [255 127 0]./255, [209 151 105]./255, ... % orange, tan
%     [251 180 185]./255, [251 180 185]./255, [251 180 185]./255, ...% soft pink
%     [247 104 161]./255,... % medium pink
%     [197 27 138]./255, [197 27 138]./255, [197 27 138]./255, ...% dark magenta
%     [122 1 119]./255,... % dark purple
%     [161 218 180]./255, [44 127 184]./255, ... % teal, medium blue
%     [65 182 196]./255, [37 52 148]./255, ...%light blue, dark blue
%     [0.3 0.3 0.3] ...   % dark gray
%     };
