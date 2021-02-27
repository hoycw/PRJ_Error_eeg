%% SEREEGA FRN N2 and RewP P3 simulation
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([app_dir 'new_eeglab/eeglab/']);
eeglab;
addpath(genpath([app_dir 'SEREEGA/']));
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Settings
trl_len  = 800;     % post-feedback length in ms
bsln_len = 200;     % pre-feedback length in ms
n_epochs = 100;     % the number of epochs to simulate
srate    = 250;     % smapling rate in Hz (250 to match EEG data)

debug = 0;

%% Generate random data
epochs = struct();
epochs.n = n_epochs;             
epochs.srate = srate;        % their sampling rate in Hz
epochs.length = bsln_len + trl_len;       % their length in ms

%   Adds time stamps, channel names, etc.
epochs.marker = 'fb_onset';  % the epochs' time-locking event marker
epochs.prestim = bsln_len;       % pre-stimulus period in ms. note that this
                            % only affects the time indicated in the final
                            % dataset; it is ignored during simulation.
                            % i.e., a simulated latency of 500 ms becomes a
                            % latency of 300 ms when a prestimulus period
                            % of 200 ms is later applied. you can use
                            % utl_shift_latency to shift all latencies in a
                            % class to fit the pre-stimulus period.

%% load headmodel 
% Add NY Head (ICBM-NY) lead field
addpath([root_dir 'PRJ_Error_eeg/data/lead_field_models/']);

% Load lead field
leadfield = lf_generate_fromnyhead('montage', 'S64');

if debug
    plot_headmodel(leadfield);
end

%% noise
noise_brown = struct( ...
    'type', 'noise', ...
    'color', 'brown', ...
    'amplitude', 2);
noise_brown = utl_check_class(noise_brown);

noise_white = struct( ...
    'type', 'noise', ...
    'color', 'white', ...
    'amplitude', 2);
noise_white = utl_check_class(noise_white);

%% Draw P3 from template
p3a_amp    = 20;
p3a_pk_lat = 380;

p3b_amp    = [14 4];
p3b_pk_lat = [300 400];

% p3_comp = utl_get_component_fromtemplate('p300_erp', leadfield);
p3a_comp = utl_get_component_fromtemplate('p3a_erp', leadfield);
p3b_comp = utl_get_component_fromtemplate('p3b_erp', leadfield);

% Modify template
p3a_comp.signal{1}.peakAmplitude = p3a_amp;
p3a_comp.signal{1}.peakLatency   = p3a_pk_lat + bsln_len;
p3b_comp.signal{1}.peakAmplitude = p3b_amp;
p3b_comp.signal{1}.peakLatency   = p3b_pk_lat + bsln_len;

plot_signal_fromclass(p3a_comp.signal{1}, epochs);
plot_signal_fromclass(p3b_comp.signal{1}, epochs);

% % Create [chan, samples, epochs] matrix
% % p3_scalpdata = generate_scalpdata(p3_comp, leadfield, epochs);
% p3a_scalpdata = generate_scalpdata(p3a_comp, leadfield, epochs);
% p3b_scalpdata = generate_scalpdata(p3b_comp, leadfield, epochs);
% 
% % Generate data in EEGlab format
% % p3_eeglab = utl_create_eeglabdataset(p3_scalpdata, epochs, leadfield);
% p3a_eeglab = utl_create_eeglabdataset(p3a_scalpdata, epochs, leadfield);
% p3b_eeglab = utl_create_eeglabdataset(p3b_scalpdata, epochs, leadfield);
% 
% % Convert to Fieldtrip
% % p3_ft = eeglab2fieldtrip(p3_eeglab, 'raw', 'none');
% p3a_ft = eeglab2fieldtrip(p3a_eeglab, 'raw', 'none');
% p3b_ft = eeglab2fieldtrip(p3b_eeglab, 'raw', 'none');

% Plot data
if debug
    cfg_plot = [];
    cfg_plot.viewmode = 'vertical';
    cfg_plot.channel = {'Fz','FCz','Cz','CPz','Pz','Oz'};
    ft_databrowser(cfg_plot, p3a_ft);
    
    ft_databrowser(cfg_plot, p3b_ft);
end

%% Theta N2
% Theta parameters
theta_freq = [3 5 6 8];
theta_amp  = 20;
theta_onset = 200;  % center of the burst in ms after feedback
theta_width = 250;  % width (half duration) of burst in ms
theta_taper = 0.5;  % taper of the burst
theta_loc_r = [11 24 24];   % coordinates from Hauser et al. 2014 NeuroImage
theta_loc_l = [-11 24 24];   % coordinates from Hauser et al. 2014 NeuroImage

% Create theta burst
theta_ersp = struct('type', 'ersp');
theta_ersp.frequency  = theta_freq;
theta_ersp.amplitude  = theta_amp;
theta_ersp.modulation = 'burst';
theta_ersp.modLatency = theta_onset + theta_width/2 + bsln_len;
theta_ersp.modWidth   = theta_width;
theta_ersp.modTaper   = theta_taper;

theta_ersp = utl_check_class(theta_ersp);

if debug
    plot_signal_fromclass(theta_ersp, epochs);
end

% Create theta source
theta_source_r = lf_get_source_nearest(leadfield, theta_loc_r);
theta_source_l = lf_get_source_nearest(leadfield, theta_loc_l);
% plot_source_location(theta_source_r, leadfield, 'mode', '3d');
% plot_source_location(theta_source_l, leadfield, 'mode', '3d');
% plot_source_projection(theta_source_r, leadfield);
% plot_source_projection(theta_source_l, leadfield);

% Compile theta component
theta_comp = struct();
theta_comp.signal = {theta_ersp, theta_ersp};
theta_comp.source  = [theta_source_r; theta_source_l];
theta_comp = utl_check_component(theta_comp, leadfield);

% % Create [chan, samples, epochs] matrix
% theta_scalpdata = generate_scalpdata(theta_comp, leadfield, epochs);
% 
% % Generate data in EEGlab format
% theta_eeglab = utl_create_eeglabdataset(theta_scalpdata, epochs, leadfield);
% 
% % Convert to Fieldtrip
% theta_ft = eeglab2fieldtrip(theta_eeglab, 'raw', 'none');

% Plot data
if debug
    cfg_plot = [];
    cfg_plot.viewmode = 'vertical';
    cfg_plot.channel = {'Fz','FCz','Cz','CPz','Pz','Oz'};
    ft_databrowser(cfg_plot, theta_ft);
end

%% Generate 25 random noise components
sources = lf_get_source_spaced(lf, 10, 25);
noise_brown.amplitude = 5;
noise_comps = utl_create_component(sources, noise_brown, lf);
noise_comps = utl_add_signal_tocomponent(noise_white,noise_comps);

% Split into 2 conditions
[comps1, comps2] = deal(noise_comps);

%% Combine Theta, P3a, P3b, and noise
comps1(end+1) = p3a_comp;
comps1(end+1) = p3b_comp;
comps1(end+1) = theta_comp;

data1 = generate_scalpdata(comps1,leadfield, epochs);
data1_eeglab = utl_create_eeglabdataset(theta_scalpdata, epochs, leadfield);
data1_ft = eeglab2fieldtrip(data1_eeglab, 'raw', 'none');

% Plot data
if debug
    cfg_plot = [];
    cfg_plot.viewmode = 'vertical';
    cfg_plot.channel = {'Fz','FCz','Cz','CPz','Pz','Oz'};
    ft_databrowser(cfg_plot, data1_ft);
end
