function [fwhm] = fn_compute_wavelet_fwhm(time_vec,freqs,n_cycles,varargin)
%% Compute the full-width at half-maximum in sec for frequencies of a wavelet signal
%   Code based on Supplement from MX Cohen 2019 NeuroImage
% INPUTS:
%   time_vec [array] - time vector for filtered data
%       e.g. [-pad_len:1/sample_rate:pad_len] where pad_len*2 is the TFR window
%   freqs [array] - vector of frequencies
%   n_cycles [array or int] - number of cycles per wavelet
% OUTPUTS:
%   fwhm [array] - full-width at half-max in sec for each frequency
%   if plot_fig, plots the gaussian windows and scatter of fwhm vs. freqs

%% Convert inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'plot_fig')
            plot_fig = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('plot_fig','var'); plot_fig = false; end

if abs(time_vec(1))~=abs(time_vec(end)); error('Time vector not centered on zero!'); end
if numel(n_cycles)~=numel(freqs) && numel(n_cycles)==1
    n_cycles = repmat(n_cycles,size(freqs));
end


% Find center of the window
win_center = dsearchn(time_vec',0);

%% Compute FWHM
fwhm  = zeros(size(freqs));
std_freq = zeros(size(freqs));
std_time = zeros(size(freqs));
gauss_std_fwhm_factor = 2.35482004503;  % ratio to convert gaussian stdev to FWHM
for f_ix = 1:length(freqs)
%     % Create the Gaussian using the FWHM formula (equation 3)
%     gauss_win = exp( (-4*log(2)*time_vec.^2) ./ fwhm(fi)^2 );
%     % measure the empirical fwhm
%     fwhm_ms(fi,1) = time_vec(midp-1+dsearchn(gauss_win(midp:end)',.5)) - time_vec(dsearchn(gauss_win(1:midp)',.5));
    
    % Fieldtrip computations:
    std_freq(f_ix) = freqs(f_ix)/n_cycles(f_ix);
    std_time(f_ix) = 1/(2*pi*std_freq(f_ix));
    ft_fwhm(f_ix) = std_time(f_ix)*gauss_std_fwhm_factor;
    
    % Create the Gaussian using the n-cycles formula (equations 1-2)
    gauss_width = n_cycles(f_ix) / (2*pi*freqs(f_ix));
    gauss_win   = exp( -time_vec.^2 ./ (2*gauss_width^2) );
    
    % empirical FWHM
    fwhm_start = time_vec(dsearchn(gauss_win(1:win_center)',.5));
    fwhm_end   = time_vec(win_center-1+dsearchn(gauss_win(win_center:end)',.5));
    fwhm(f_ix) = fwhm_end - fwhm_start;
end

%% Plot Results
if plot_fig
    figure; hold on;
    f_lines  = gobjects(size(freqs));
    f_leg    = cell(size(freqs));
    cmap_idx = zeros(size(freqs));
    for f_ix = 1:numel(freqs)
        % Plot Gaussian window
        subplot(1,2,1); hold on;
        cmap_idx(f_ix) = f_ix*floor(size(cmap,1)/numel(freqs));%dsearchn([1:size(cmap,1)]'./numel(freqs),f_ix);
        f_lines(f_ix) = plot(time_vec,gauss_win,'Color',cmap(cmap_idx(f_ix),:));
        f_leg{f_ix} = [num2str(freqs(f_ix)) ' Hz: ' num2str(fwhm(f_ix))];
    end
    legend(f_lines,f_leg);
    
    % Plot FWHM metrics
    subplot(1,2,2); hold on;
    mxc = scatter(freqs,fwhm,50,'r');
    ft = scatter(freqs,ft_fwhm,50,'b');
    legend([mxc mxc_h ft], {'Cohen', 'Cohen H', 'FT'});
end

%%
% figure(3), clf
% 
% plot(freqs,fwhm_ms*1000,'o-','markersize',8,'markerfacecolor','w','linew',2)
% xlabel('Wavelet frequency (Hz)'), ylabel('FWHM (ms)')
% legend({'Using FWHM';'Using n-cycles'})
% set(gca,'xlim',[freqs(1)-1 freqs(end)+1])


end