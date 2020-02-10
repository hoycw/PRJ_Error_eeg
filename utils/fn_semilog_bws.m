function [bws] = fn_semilog_bws(freqs)
%% Compute the bandwidth of a Gaussian filter bank based on semilog spacing
% From Erik Edward's thesis:
%   in log(BWstd) = 0.5*log(f) + 0.39;
% INPUTS:
%   freqs [float] - array of frequencies for which adesired filter bandwidths
% OUTPUTS:
%   bws [float] - bandwidths for each frequency, returned as standard
%       deviation of the gaussian filter in frequency space

bws = zeros(size(freqs));
for f = 1:numel(freqs)
    bws(f) = pi^(0.5*log(freqs(f)) + 0.39);
end

end