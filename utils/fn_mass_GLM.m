function betas = fn_mass_GLM(X,y,add_offset)
%function betas = massGLM(X,y,add_offset)
%efficiently calculates betas for mass univariate GLM
%
%inputs: X (design matrix) should be nTrials*nRegressors
%        y (data) should be nTrials*nChannels*nTimepoints(*nFrequencies)
%        add_offset should be true to add a constant regressor to
%        the design matrix, false if not
%
%outputs: betas is nRegressors*nChannels*nTimepoints(*nFrequencies)

if nargin < 3
    add_offset = true;
end

%make sure inputs are in order
if size(X,2)>size(X,1)
    X = X';
    warning('massGLM: transposing design matrix')
end
if any(all(X==1,1)) & add_offset
    add_offset = false;
    warning('massGLM: design matrix already appears to have a constant regressor.')
end

%get data and design matrix dimensions
dims  = size(y);
ntrl  = size(X,1);
if dims(1) ~= ntrl
    error('massGLM: trial numbers in X and y don''t match!');
end
X     = [ones(ntrl,add_offset) X];
nregs = size(X,2);

%reshape data matrix so it is nTrials*nFeatures
%(=channelsXtimepointsXfrequencies...)
y  = reshape(y,[dims(1) prod(dims(2:end))]);

%estimate betas
pX    = pinv(X);
betas = pX*y;

%reshape betas to sensible output
betas = reshape(betas,[nregs dims(2:end)]);
end