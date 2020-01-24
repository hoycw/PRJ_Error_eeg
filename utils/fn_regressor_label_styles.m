function [reg_labels, colors, line_styles] = fn_regressor_label_styles(model_id)
%% Converts the name of a model into regressor labels, plotting colors/styles
% Future:
%   add unsigned distance metrics
%   add categorical 0/1 if becomes necessary
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

%% List of possible regressors and their colors
regressors  = {'pWin','sPE','uPE','tRT','ptRT','p2tRT','lDist','ulDist','iDist'};
reg_colors = {[152 78 163]./255, [255 127 0]./255, [166 86 40]./255,...
    [247 129 191]./255, [247 129 191]./255, [247 129 191]./255, ...
    [0.5 0.5 0.5], [0.5 0.5 0.5], [0.5 0.5 0.5]};
% Color options:
%   purple: [152 78 163]
%   orange: [255 127 0]
%   brown:  [166 86 40]
%   yellow: [255 255 51]
%   pink:   [247 129 191]
%   red:    [228 26 28]
%   green:  [55 126 184]
%   blue:   [77 175 74]
%   gray:   [0.5 0.5 0.5]
% Taken from: colorbrewer2.org, qualitative, 5-class Set1


%% Convert model_id into set of conditions
switch model_id
    case {'RL','pWinPEus'}
        reg_labels = {'pWin','sPE','uPE'};
    case 'RLRT'
        reg_labels = {'pWin','sPE','uPE','tRT'};
    case 'RLpRT'
        reg_labels = {'pWin','sPE','uPE','tRT','ptRT'};
    case 'RLpRTlD'
        reg_labels = {'pWin','sPE','uPE','tRT','ptRT','lDist'};
    case 'RLpRTulD'
        reg_labels = {'pWin','sPE','uPE','tRT','ptRT','lDist','ulDist'};
    case 'RLpRTiD'
        reg_labels = {'pWin','sPE','uPE','tRT','ptRT','iDist'};
    otherwise
        error(strcat('Unknown model_id: ',model_id));
end

% Assign colors and line styles
colors = cell(size(reg_labels));
line_styles = cell(size(reg_labels));
for reg_ix = 1:numel(reg_labels)
    colors{reg_ix} = reg_colors{strcmp(reg_labels{reg_ix},regressors)};
    if any(strcmp(reg_labels{reg_ix},{'ptRT','ulDist'}))
        line_styles{reg_ix} = ':';
    elseif strcmp(reg_labels{reg_ix},'p2tRT')
        line_styles{reg_ix} = '-.';
    else
        line_styles{reg_ix} = '-';
    end
end

end
