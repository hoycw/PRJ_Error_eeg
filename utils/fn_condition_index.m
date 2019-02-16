function condition_n = fn_condition_index(conditions, bhv)
% Returns index of trial condition assignments based on requested conditions
% INPUTS:
%   conditions [str] - name of the set of conditions requested
%   bhv [struct] - trial info structure containing info for logic sorting
%       Total_Trial,Block,Condition,Hit,RT,Timestamp,Tolerance,Trial,Score,ITI,ITI type
% OUTPUTS:
%   condition_n [int vector] - integer assignment of each trial based on conditions

[cond_lab, ~, ~, ~] = fn_condition_label_styles(conditions);

condition_n = zeros(size(bhv.Trial));
for cond_ix = 1:numel(cond_lab)
    switch cond_lab{cond_ix}
        case 'Ez'
            condition_n(strcmp('easy',bhv.Condition)) = cond_ix;
        case 'Hd'
            condition_n(strcmp('hard',bhv.Condition)) = cond_ix;
        case 'Wn'
            condition_n(bhv.Hit==1) = cond_ix;
        case 'Ls'
            condition_n(bhv.Hit==0) = cond_ix;
        case 'Er'
            warning('WARNING!!! Assuming target_time = 1 sec...');
            condition_n(bhv.RT<1) = cond_ix;
        case 'Lt'
            warning('WARNING!!! Assuming target_time = 1 sec...');
            condition_n(bhv.RT>1) = cond_ix;
        case 'EzWn'
            matches = logical(strcmp('easy',bhv.Condition)) & logical(bhv.Hit==1);
            condition_n(matches) = cond_ix;
        case 'EzLs'
            matches = logical(strcmp('easy',bhv.Condition)) & logical(bhv.Hit==0);
            condition_n(matches) = cond_ix;
        case 'HdWn'
            matches = logical(strcmp('hard',bhv.Condition)) & logical(bhv.Hit==1);
            condition_n(matches) = cond_ix;
        case 'HdLs'
            matches = logical(strcmp('hard',bhv.Condition)) & logical(bhv.Hit==0);
            condition_n(matches) = cond_ix;
        otherwise
            error(['Invalid condition label: ' conditions]);
    end
end

if sum(condition_n==0)~=0
    warning(['Not all trials accounted for by conditions: ' conditions]);
end

end
