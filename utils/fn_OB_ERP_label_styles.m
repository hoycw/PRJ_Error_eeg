function [labels, names, colors, line_styles, markers] = fn_OB_ERP_label_styles(grp_id)
%% Converts the name of a group of OB ERP features into labels, plotting colors/styles
% INPUTS:
%   grp_id [str] - label for group of conditions to select
%       'N2b': N200 from Novelty 'Odd'
%       'P3a': P300 from Novelty 'Odd'
%       'N2c': N200 from Target 'Tar'
%       'P3b': P300 from Target 'Tar'
% OUTPUTS:
%   labels [cell array] - string short-hand labels of specific conditions
%   names [cell array] - string longer, full labels of specific conditions
%   colors [cell array] - [R G B] tuples (scaled to 0-1) per condition
%   line_styles [cell array] - solid ('-') or dotted ('--') per condition
%   markers [cell array] - {'*','o','d'} markers per condition

if ~strcmp(grp_id,'N2sP3s'); error('N2sP3s is only valid model right now!'); end

%% List of possible labels and their colors
% Condition codes
features  = {'N2b','P3a','N2c','P3b'};

% Condition labels (longer strings for clear plotting legends)
condition_names = {'Novel N2b','Novel P3a','Target N2c','Target P3b'};

% Condition Colors ([R G B] scaled to 0-1)
condition_colors = {...
    [197,27,138]./256, [197,27,138]./256, ...      % Novel Magenta
    [56 108 176]./256, [56 108 176]./256 ...       % Target Cyan
%     [158 188 218]./256;...  % lighter purple
%     [254 204 92]./256;...   % lighter red
%     [129 15 124]./256;...   % darker purple
%     [189 0 38]./256 ...     % darker red
    };

%% Convert grp_id into set of conditions
switch grp_id
    % ---------------------------------------------------------------------
    % Oddball ERPs
    case 'N2sP3s'
        labels = {'N2b','P3a','N2c','P3b'};
    case 'N2s'
        labels = {'N2b','N2c'};
    case 'P3s'
        labels = {'P3a','P3b'};
    
    otherwise
        error(strcat('Only one, unrecognized condition offered: ',grp_id));
end

%% Assign colors and line styles
names       = cell(size(labels));
colors      = cell(size(labels));
line_styles = cell(size(labels));
markers     = cell(size(labels));
for cond_ix = 1:numel(labels)
    names{cond_ix}  = condition_names{strcmp(labels{cond_ix},features)};
    colors{cond_ix} = condition_colors{strcmp(labels{cond_ix},features)};
    
    % Define Line Styles
    if any(strcmp(labels{cond_ix},{'N2b','N2c'}))
        line_styles{cond_ix} = '-';
    elseif any(strcmp(labels{cond_ix},{'P3a','P3b'}))
        line_styles{cond_ix} = '--';
    else
        error(['What OB ERP feature is ' labels{cond_ix} '?']);
    end
    
    % Define Marker Styles
    if any(strcmp(labels{cond_ix},{'N2b','N2c'}))
        markers{cond_ix} = '*';     % Surprise (difficulty agnostic)
    elseif any(strcmp(labels{cond_ix},{'P3a','P3b'}))
        markers{cond_ix} = 'o';     % Easy
    end
end

end
