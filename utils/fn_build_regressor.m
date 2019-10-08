function [regressor] = fn_build_regressor(factor, bhv)%contrast
%% Build regressors for GLM analysis
% INPUTS:
%   factor [str] - name of the regressor
%   bhv [struct] - behavioral data to find trial matches
%   FUTURE:
%       contrast [0/1] - make the regressor a contrast (-1/1) or not (binary 0/1)
% OUTPUTS:
%   regressor [int matrix] - contrast (-1/1) or binary (0/1) regressor of interest

regressor = nan([numel(bhv.trl_n) 1]);
if strcmp(factor,'Out')
    regressor(fn_condition_index({'Wn'}, bhv)==1) = 1;
    regressor(fn_condition_index({'Ls'}, bhv)==1) = -1;
elseif strcmp(factor,'Dif')
    regressor(fn_condition_index({'Hd'}, bhv)==1) = 1;
    regressor(fn_condition_index({'Ez'}, bhv)==1) = -1;
elseif strcmp(factor,'Sur')
    regressor(fn_condition_index({'EzLs','HdWn','Su'}, bhv)~=0) = 1;
    regressor(fn_condition_index({'EzWn','HdLs'}, bhv)~=0) = -1;
elseif strcmp(factor,'off')
    regressor = ones(size(bhv.trl_n));
else
    error(['Unknown factor: ' factor]);
end

if any(isnan(regressor))
    fprintf(2,'%i / %i unassigned trials in factor %s! Assigning 0...\n',sum(isnan(regressor)),numel(bhv.trl_n),factor);
    regressor(isnan(regressor)) = 0;
end

end