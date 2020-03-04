function condition_n = fn_condition_index(cond_lab, bhv)
% Returns index of trial condition assignments based on requested conditions
% INPUTS:
%   cond_lab [cell array] - strings with names of the set of conditions requested
%   bhv [struct] - trial info structure containing info for logic sorting
%       Total_Trial, Block, Feedback, RT, Timestamp, Tolerance, Trial, Hit, Score, bad_fb, Condition, ITI, ITI type
% OUTPUTS:
%   condition_n [int vector] - integer assignment of each trial based on conditions

condition_n = zeros(size(bhv.trl_n));
for cond_ix = 1:numel(cond_lab)
    if any(strcmp(cond_lab{cond_ix},{'ErQ1','ErQ2','LtQ3','LtQ4'}))
        quartiles = quantile(bhv.rt,4);
    end
    switch cond_lab{cond_ix}
        case 'Ez'
            condition_n(strcmp('easy',bhv.cond)) = cond_ix;
        case 'Hd'
            condition_n(strcmp('hard',bhv.cond)) = cond_ix;
        case 'Wn'
            condition_n(strcmp(bhv.fb,'W')) = cond_ix;
        case 'Ls'
            condition_n(strcmp(bhv.fb,'L')) = cond_ix;
        case 'Su'
            condition_n(strcmp(bhv.fb,'S')) = cond_ix;
        case 'Er'
            warning('WARNING!!! Assuming target_time = 1 sec...');
            condition_n(bhv.rt<1) = cond_ix;
        case 'Lt'
            warning('WARNING!!! Assuming target_time = 1 sec...');
            condition_n(bhv.rt>1) = cond_ix;
        case 'ErQ1'
            condition_n(bhv.rt<=quartiles(1)) = cond_ix;
        case 'ErQ2'
            condition_n(bhv.rt>quartiles(1) &bhv.rt<=quartiles(2)) = cond_ix;
        case 'LtQ3'
            condition_n(bhv.rt>quartiles(2) & bhv.rt<=quartiles(3)) = cond_ix;
        case 'LtQ4'
            condition_n(bhv.rt>quartiles(3)) = cond_ix;
        case 'EzWn'
            matches = logical(strcmp('easy',bhv.cond)) & strcmp(bhv.fb,'W');
            condition_n(matches) = cond_ix;
        case 'EzLs'
            matches = logical(strcmp('easy',bhv.cond)) & strcmp(bhv.fb,'L');
            condition_n(matches) = cond_ix;
        case 'EzSu'
            matches = logical(strcmp('easy',bhv.cond)) & strcmp(bhv.fb,'S');
            condition_n(matches) = cond_ix;
        case 'HdWn'
            matches = logical(strcmp('hard',bhv.cond)) & strcmp(bhv.fb,'W');
            condition_n(matches) = cond_ix;
        case 'HdLs'
            matches = logical(strcmp('hard',bhv.cond)) & strcmp(bhv.fb,'L');
            condition_n(matches) = cond_ix;
        case 'HdSu'
            matches = logical(strcmp('hard',bhv.cond)) & strcmp(bhv.fb,'S');
            condition_n(matches) = cond_ix;
        case 'Odd'
            matches = logical(strcmp('odd',bhv.cond));
            condition_n(matches) = cond_ix;
        case 'Std'
            matches = logical(strcmp('std',bhv.cond));
            condition_n(matches) = cond_ix;
        case 'Tar'
            matches = logical(strcmp('tar',bhv.cond));
            condition_n(matches) = cond_ix;
        otherwise
            error(['Invalid condition label: ' cond_lab]);
    end
end

if sum(condition_n==0)~=0
    warning(['Not all trials accounted for by conditions: ' strjoin(cond_lab,',')]);
end

end
