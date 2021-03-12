function [labels, names, colors, line_styles, markers] = fn_regressor_label_styles(model_id)
%% Converts the name of a model into regressor labels, plotting colors/styles
%   ^ = regressors were only used for early model comparison (not in the paper)
%       Generally, RL features are better than anything else for nearly all analyses
%   RL Features: expected value (logistic regression), signed/unsigned RPE
%   Outcome probability/likelihood
%   ^Categorical Reward/Outcome Features: Valence/Sign, Value, Magnitude
%   ^Expected Value alternatives: block accuracy, rolling accuracy (5 or 10 trials)
%   design features: inter-trial interval (ITI)
%   statistical controls: SBJ (random intercepts in LME)
% INPUTS:
%   model_id [str] - format is ['model label' '_' 'cond_lab']
%       Final paper model: 'ERPEsL_DifFB' ~ EV + sRPE + uRPE + Lik + (1|SBJ)
% OUTPUTS:
%   labels [cell array] - string short-hand labels of specific regressors
%   names [cell array] - string longer, full labels of specific regressors
%   colors [cell array] - [R G B] tuples (scaled to 0-1) per regressor
%   line_styles [cell array] - line plotting styles per regressor
%   markers [cell array] - scatter plot markers per regressor
% Color options: Taken from colorbrewer2.org, qualitative, 5-class Set1
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

%% List of possible regressors and their colors
regressors  = {...
    'Sign','Val','Mag',... % Categorical reward features
    'EV','bAcc','rAcc','rAcc10','score',... % Outcome predictors
    'sRPE','uRPE',... % Reward Feedback
    'Lik','rLik',... % Outcome Likelihood, residualized Likelihood
    'ERB','rough','mxS','mnS','mxdS','mndS', ... % Auditory salience properties (via Sijia Zhao)
    'ITI', ... % preceding inter-trial interval
    'SBJ' ... % SBJ only null model (baseline model fit using only SBJ random intercepts)
    };
regressor_names = {...
    'Valence','Value','Magnitude',...
    'Expected Value','Block Accuracy','Rolling Accuracy','Rolling Accuracy 10','Total Score',...
    'RPE Value','RPE Magnitude',...
    'Likelihood','Likelihood',...
    'Loudness','Roughness','Max Salience','Mean Salience','Max dSalience','Mean dSalience',...
    'ITI', ...
    'SBJ' ...
    };
regressor_colors = {...
    [209 151 105]./255, [209 151 105]./255, [118 160 156]./255, ... % tan, tan, teal, purple [152 78 163]./255,
    [0 0 0], [0 0 0], [0 0 0], [0 0 0], [0.4 0.4 0.4],... % dark blue, blacks, gray
    [209 151 105]./255, [118 160 156]./255, ... % tan, teal
    [152 78 163]./255,[152 78 163]./255,... % purple
    [231,41,138]./255, [102,166,30]./255, ... %magenta, lime green
    [217,95,2]./255, [217,95,2]./255, [117,112,179]./255, [117,112,179]./255, ... % dark orange, mauve purple
    [0.3 0.3 0.3], ...   % dark gray
    [0.3 0.3 0.3] ...   % dark gray
    };
regressor_markers = {...
    's','s','d',...%'*',
    'x','x','x','x','^',...
    's','d',...
    '*','*',...
    '*','*','*','*','*','*',...
    'v','v'...
    };

%% Convert model_id into set of conditions
switch model_id
    %======================================================================
    % Categorical Outcome Features
    case 'S'
        labels = {'Sign'};
    case 'V'
        labels = {'Val'};
    case 'VL'
        labels = {'Val','Lik'};
    case 'VM'
        labels = {'Val','Mag'};
    case 'ML'
        labels = {'Mag','Lik'};
    case 'VML'
        labels = {'Val','Mag','Lik'};
    case 'SML'
        labels = {'Sign','Mag','Lik'};
    case 'VSML'
        labels = {'Val','Sign','Mag','Lik'};
    
    %======================================================================
    % RL models without EV
    case 'sRPE'
        labels = {'sRPE'};
    case 'uRPE'
        labels = {'uRPE'};
    case 'RPEs'
        labels = {'sRPE','uRPE'};
    case 'RPEsL'
        labels = {'sRPE','uRPE','Lik'};
    case 'uRPEL'
        labels = {'uRPE','Lik'};
        
    case 'RSVPE'
        labels = {'Sign','Val','sRPE'};
        
    %======================================================================
    % RL models with EV
    case 'ERPEs'
        labels = {'EV','sRPE','uRPE'};
    case 'EsRPE'
        labels = {'EV','sRPE'};
    case 'EsRPEL'
        labels = {'EV','sRPE','Lik'};
    case 'EuRPE'
        labels = {'EV','uRPE'};
    case 'ERPEsL'
        labels = {'EV','sRPE','uRPE','Lik'};
    case 'ERPEsrL'
        labels = {'EV','sRPE','uRPE','rLik'};
    case 'ERPEsiti'
        labels = {'EV','sRPE','uRPE','ITI'};
    case 'ERPEsscr'
        labels = {'EV','score','sRPE','uRPE'};
        
    %======================================================================
    % Combinations of Outcome and RL models
    case 'VMLsRPE'
        labels = {'Val','Lik','Mag','sRPE'};
    case 'VMLRPEsL'
        labels = {'Val','Lik','Mag','EV','sRPE','uRPE','Lik'};
        
    %======================================================================
    % Auditory Salience Features
    %   Only ERB should be used because "loudness" is different across my
    %   sounds, making the other features uninterpretable (even if they fit)
    case 'ERB'
        labels = {'ERB'};
%     case 'rough'
%         labels = {'rough'};
%     case 'ERBr'
%         labels = {'ERB','rough'};
%     case 'Kayser'
%         labels = {'mxS','mnS','mxdS','mndS'};
%     case 'AudSal'
%         labels = {'ERB','rough','mxS','mnS','mxdS','mndS'};
    case 'ERBuRPE'
        labels = {'uRPE','ERB'};
    case 'ERBsRPE'
        labels = {'sRPE','ERB'};
%     case 'rsRPE'
%         labels = {'sRPE','rough'};
    case 'ERBrsRPE'
        labels = {'sRPE','ERB','rough'};
%     case 'ASsRPE'
%         labels = {'sRPE','ERB','rough','mxS','mnS','mxdS','mndS'};
    
    %======================================================================
    % RL models with alternative accuracy
%     case 'RLbA'
%         labels = {'bAcc','sRPE','uRPE'};
%     case 'RLbApW'
%         labels = {'bAcc','pWin','sRPE','uRPE'};
%     case 'RLrA10'
%         labels = {'rAcc10','sRPE','uRPE'};
    case 'rARPEs'
        labels = {'rAcc','sRPE','uRPE'};
%     case 'RLrApW'
%         labels = {'rAcc','pWin','sRPE','uRPE'};
    
    %======================================================================
    % Null SBJ only control model
    case 'SBJonly'
        labels = {'SBJ'};
        
    %======================================================================
    % Previous Trial Regressors:
%     case 'RLpT'
%         error('why running with previous trial regressors?');
%         labels = {'pWin','sRPE','uRPE','sTaEr','psTaEr'};
%     case 'RLpRTlD'
%         error('why running with previous trial regressors?');
%         labels = {'pWin','sRPE','uRPE','sTaEr','psTaEr','sThPr'};
%     case 'RLpRTulD'
%         error('why running with previous trial regressors?');
%         labels = {'pWin','sRPE','uRPE','sTaEr','psTaEr','sThPr','uThPr'};
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
markers     = cell(size(labels));
for reg_ix = 1:numel(labels)
    names{reg_ix}   = regressor_names{strcmp(labels{reg_ix},regressors)};
    colors{reg_ix}  = regressor_colors{strcmp(labels{reg_ix},regressors)};
    markers{reg_ix} = regressor_markers{strcmp(labels{reg_ix},regressors)};
    line_styles{reg_ix} = '-';
%     if any(strcmp(labels{reg_ix},'psTaEr'))%{'psTaEr','uThPr'} % strcmp(labels{reg_ix}(1),'u') || 
%         line_styles{reg_ix} = ':';
%     elseif strcmp(labels{reg_ix},'p2sTaEr')
%         line_styles{reg_ix} = '-.';
%     else
%         line_styles{reg_ix} = '-';
%     end
end

% Check for overlap in type of regressor
unique_colors = unique(vertcat(colors{:}),'rows');
if size(unique_colors,1)~=numel(labels)
    for reg_ix = 1:numel(labels)
        if any(strcmp(labels{reg_ix},{'Val','Mag','mnS','mndS'}))
            line_styles{reg_ix} = ':';
        elseif strcmp(labels{reg_ix},'Sign')
            line_styles{reg_ix} = '-.';
        end
    end
end

end
