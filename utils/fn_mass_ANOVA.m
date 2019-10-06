function omega = fn_mass_ANOVA(y,group)
%function omega = fn_mass_ANOVA(y,group)
%efficiently calculates omega squared values for mass univariate ANOVA
%
%inputs: group (design matrix) should be a cell array of nTrials*1 factors
%        y (data) should be nTrials*nChannels*nTimepoints(*nFrequencies)
%        
%outputs: omega is nFactors*nChannels*nTimepoints(*nFrequencies)

%get data and design matrix dimensions
dims  = size(y);
ntrl  = size(group{1},1);
if dims(1) ~= ntrl
    error('massANOVA: trial numbers in group and y don''t match!');
end
nfac = size(group,2); %number of factors in ANOVA

%reshape data matrix so it is nTrials*nFeatures
%(=channelsXtimepointsXfrequencies...)
y  = reshape(y,[dims(1) prod(dims(2:end))]);

%run ANOVA
table = anova2D(y,group);

%estimate omega-squared
mse = table.mse; %mean squared error
sst = table.sst; %total sum of squares
w   = nan(size(table.ssterm),'single');
for icond = 1:size(table.ssterm,1) %loop over factors
    df  = table.dfterm(icond); %d.f. for effect of interest
    ss  = table.ssterm(icond,:);%sum of squares for effect of interest
    w(icond,:) = (ss - (df .* mse)) ./ (sst + mse); %omega squared
end
        
%reshape omega squareds to sensible output dimensions
omega = reshape(w,[nfac dims(2:end)]);
end
