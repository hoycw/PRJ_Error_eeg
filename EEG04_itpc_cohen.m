% mikexcohen@gmail.com

%% Compute and plot TF-ITPC for one electrode
function itpc(filename, SBJ)
addpath(genpath('/Users/SCS22/Desktop/Knight_Lab/Preprocessing_Work'));
load(filename,'data');
for x= 1: numel(data.label);
% wavelet parameters
num_frex = 40;
min_freq =  2;
max_freq = 30;

channel = data.label(x);

% set range for variable number of wavelet cycles
range_cycles = [ 4 10 ];

% other wavelet parameters
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
wavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);
time = -2:1/data.fsample:2;
half_wave_size = (length(time)-1)/2;

% FFT parameters
nWave = length(time);
nData = length(data.time)*326;
nConv = nWave+nData-1;
reshappedfft = zeros(1,nData);
for n = 1: length(data.trial);
    for m = 1:326;
    reshapedfft(1,(n+(n-1)*m)) = data.trial{n}(x,m);
    end
end

% FFT of data (doesn't change on frequency iteration)
dataX = fft(reshapedfft,nConv);

% initialize output time-frequency data
global tf
tf = zeros(num_frex,326);

% loop over frequencies
for fi=1:num_frex;
    
    % create wavelet and get its FFT
    s = wavecycles(fi)/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    waveletX = fft(wavelet,nConv);
    
    % run convolution
    as = ifft(waveletX.*dataX,nConv);
    as = as(half_wave_size+1:end-half_wave_size);
    as = reshape(as,326,length(data.trial));
    
    % compute ITPC
    tf(fi,:) = abs(mean(exp(1i*angle(as)),2));
end
% plot results
figure(1), clf
contourf(data.time{x},frex,tf,40,'linecolor','none')
set(gca,'clim',[0 .6],'ydir','normal','xlim',[-300 1000])
title(['ITPC - ' data.label{x}])
end

cfg = [];
cfg.channels = {'all', '-EXG*'};
cfg.layout    = 'biosemi64.lay';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, tf)
savefig(data_out_filename);
save(data_out_filename_itpc, 'tf');

%% illustration of ITPC with different variances

% specify parameters
circ_prop = .5; % proportion of the circle to fill
N = numel(data.trials); % number of "trials"

% generate phase angle distribution
simdata = rand(1,N) * (2*pi) * circ_prop;


% compute ITPC and preferred phase angle
itpc      = abs(mean(exp(1i*simdata)));
prefAngle = angle(mean(exp(1i*simdata)));


% and plot...
figure(2), clf

% as linear histogram
subplot(3,3,4)
hist(simdata,20)
xlabel('Phase angle'), ylabel('Count')
set(gca,'xlim',[0 2*pi])
title([ 'Observed ITPC: ' num2str(itpc) ])

% and as polar distribution
subplot(1,2,2)
polar([zeros(1,N); simdata],[zeros(1,N); ones(1,N)],'k')
hold on
h = polar([0 prefAngle],[0 itpc],'m');
set(h,'linew',3)
title([ 'Observed ITPC: ' num2str(itpc) ])

%% end.
