function [vifs] = fn_variance_inflation_factor(design)
%% Compute Variance Inflation Factor (VIF)
%   quantifies how much the variance is inflated due to collinearity of regressor matrix columns

R0   = corrcoef(design); % correlation matrix
vifs = diag(inv(R0))';

end