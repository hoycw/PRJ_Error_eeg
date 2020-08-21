function [labels, names, colors, line_styles, markers] = fn_condition_label_styles(grp_id)
%% Converts the name of a group of conditions into labels, plotting colors/styles
%   Mostly for Target Time, but also some Oddball conditions
%   NOTE: S commonly stands for "surprise", which is now referred to as "neutral"
% INPUTS:
%   grp_id [str] - label for group of conditions to select
%       common options: {'EH'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3
%   light red: [251 154 153]
%   dark red: [227 26 28]
%   light blue: [166 206 227]
%   dark blue: [31 120 180]
%   light green: [178 223 138]
%   dark green: [51 160 44]

%% List of possible labels and their colors
% Condition codes
conditions  = {...
    'Std','Tar','Odd',...                           % Oddball conditions
    'Ez','Hd','Wn','Ls','Su',...                    % Target Time block and outcome types
    'EzWn','EzLs','EzSu','HdWn','HdLs','HdSu',...   % Target Time conditions
    'All','Pos','Neg',...                           % Combinations of Target Time conditions
    'Er','Lt','ErQ1','ErQ2','MdQ3','LtQ4','LtQ5'};  % Performance-based target time selections

% Condition labels (longer strings for clear plotting legends)
condition_names = {...
    'Standard','Target','Oddball',...
    'Easy','Hard','Win','Loss','Surprise',...
    'Easy Win','Easy Loss','Easy Neutral','Hard Win','Hard Loss','Hard Neutral',...
    'All Trials', 'Positive', 'Negative',...
    'Early','Late','Early Q1','Early Q2','Middle Q3','Late Q4','Late Q5'};

% Condition Colors ([R G B] scaled to 0-1)
condition_colors = {...
    [168 180 165]./256, [144 205 229]./256, [142 82 126]./256, ...
    [228,26,28]./256, [55,126,184]./256, [51 160 44]./256, [227 26 28]./256, [31 120 180]./256, ...
    [51 160 44]./256, [227 26 28]./256, [31 120 180]./256, ... %lighter colors: [178 223 138]./256, [251 154 153]./256, [166 206 227]./256, ...
    [51 160 44]./256, [227 26 28]./256, [31 120 180]./256, ...
    [0 0 0], [31 120 180]./256, [227 26 28]./256,...
    [166 97 26]./256, [1 133 113]./256, ... % brown, aqua
    [166 97 26]./256, [223 194 125]./256, [0.4 0.4 0.4], [128 205 193]./256, [1 133 113]./256 ...    % [brown, tan, gray, teal, aqua]
    };

%% Convert grp_id into set of conditions
switch grp_id
    % ---------------------------------------------------------------------
    % Oddball Conditions
    case 'Odd'
        labels = {'Std','Tar','Odd'};
    case 'DifOdd_OS'
        labels = {'Odd', 'Std'};
    case 'DifOdd_TS'
        labels = {'Tar','Std'};
    case 'DifOdd_OT'
        labels = {'Odd', 'Tar'};
    
    % ---------------------------------------------------------------------
    % Target Time Conditions
    case 'All'
        labels = {'All'};
    case 'Dif'                          % Difficulty
        labels = {'Ez', 'Hd'};
    case {'OutS','FB'}                  % Feedback (includes neutral)
        labels = {'Wn', 'Ls', 'Su'};
    case 'Out'                          % Outcome (no neutral)
        labels = {'Wn', 'Ls'};
    case 'Val'                          % Valence
        labels = {'Pos', 'Neg'};
    
    % Subsets of Target Time conditions
    case 'EHSu'                         % Neutral Outcomes in Easy+Hard
        labels = {'EzSu','HdSu'};
    case {'DifOut','DifOutUE','DifOutWL','DifOutdO','DifOutDO','Holroyd'}   % Four main conditions
        labels = {'EzWn', 'EzLs', 'HdWn', 'HdLs'};
    case {'DifOutS', 'DifOutSx'}        % Four main conditions + all neutral outcomes combined
        labels = {'EzWn', 'EzLs', 'HdWn', 'HdLs', 'Su'};
    case {'DifOutSur','DifFB'}          % All 6 main conditions
        labels = {'EzWn', 'EzSu', 'EzLs', 'HdWn', 'HdSu', 'HdLs'};
    case 'EzOut'                        % Win/Loss in Easy
        labels = {'EzWn', 'EzLs'};
    case 'EzOutS'                       % All outcomes in Easy
        labels = {'EzWn', 'EzLs', 'EzSu'};
    case 'HdOut'                        % Win/Loss in Hard
        labels = {'HdWn', 'HdLs'};
    case 'HdOutS'                       % All outcomes in Hard
        labels = {'HdWn', 'HdLs', 'HdSu'};
    
    % Performance (RT) based trial selection
    case 'Tar2'                         % Performance split of early/late
        labels = {'Er', 'Lt'};
    case 'Tar5'                         % RT-based split into quintiles
        labels = {'ErQ1','ErQ2','MdQ3','LtQ4','LtQ5'};
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
        line_styles{cond_ix} = '--';%'-.';%    % Hard
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
