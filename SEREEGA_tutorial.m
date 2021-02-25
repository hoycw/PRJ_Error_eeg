%% SEREEGA Tutorial
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([app_dir 'new_eeglab/eeglab/']);
eeglab;
addpath(genpath([app_dir 'SEREEGA/']));
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Generate random data
epochs = struct();
epochs.n = 100;             % the number of epochs to simulate
epochs.srate = 1000;        % their sampling rate in Hz
epochs.length = 1000;       % their length in ms

%% Load pre-computed lead field model
% Add NY Head (ICBM-NY) lead field
addpath([root_dir 'PRJ_Error_eeg/data/lead_field_models/']);

% Select Montage
%   'S64' selects 64 channel montage from 10-20 system
%   utl_get_montage has all available montages
lf = lf_generate_fromnyhead('montage', 'S64');

% Alternative with only midline:
% lf2 = lf_generate_fromnyhead('labels',{'Fz','Cz','Pz'});

% Altenrative generated from fieldtrip:
% lf2 = lf_generate_fromfieldtrip('montage','S64','resolution',5);

% Inspect channel locations in lead field in relation to brain sources by plotting:
plot_headmodel(lf);
plot_headmodel(lf, 'style', 'boundary', 'labels', 0);

%% Pick source location
% Random source
% source = lf_get_source_random(lf);
% plot_source_location(source, lf);

% Right Visual cortex
source = lf_get_source_nearest(lf, [20 -75 0]);
plot_source_location(source, lf, 'mode', '3d');

%% Orieting source dipole
% X, Y, Z orientations
plot_source_projection(source, lf, 'orientation', [1 0 0], 'orientedonly', 1);
plot_source_projection(source, lf, 'orientation', [0 1 0], 'orientedonly', 1);
plot_source_projection(source, lf, 'orientation', [0 0 1], 'orientedonly', 1);

% NY head model default is perpendicular to cortical surface
plot_source_projection(source, lf, 'orientation', [1, 1, 0]);
plot_source_projection(source, lf);

% Fieldtrip and Pediatric have no meaningful default orientation, so can be set using:
%	utl_get_orientation_random - random orientation
%	utl_get_orientation_pseudoperpendicular - points outward toward scalp surface
%	utl_get_orientation_pseudotangential- attempted approximate tanget to cortical surface

% To match literature, must match topo, so either:
%   1) pick source with default orientation that reproduces the topo
%   2) pick a source and adjust orientation to match topo

%% Define source activation signal
% ERPs defined by latency, width, and amplitude of it's peak
erp = struct();
erp.peakLatency = 500;      % in ms, starting at the start of the epoch
erp.peakWidth = 200;        % in ms
erp.peakAmplitude = 1;      % in microvolt

% amplitude is at source, so scalp amplitude depends on lead field and
%   interference from other signals

% Check everything necessary is set
erp = utl_check_class(erp, 'type', 'erp');

% Plot ERP
plot_signal_fromclass(erp, epochs);

% Example with 2 peaks:
% erp = struct('peakLatency', [450, 500], 'peakWidth', [200, 200], 'peakAmplitude', [-1, 1]);

%% Combine ERP signal and source location into brain component
c = struct();
c.source = source;      % obtained from the lead field, as above
c.signal = {erp};       % ERP class, defined above

c = utl_check_component(c, lf);

% Can use utl_create_component for library of pre-made signal/source combos

% Options to change orientation:
% c.orientation = [0 1 0];
% c.orientation = utl_get_orientation_pseudoperpendicular(source, lf);

%% Simulate scalp data
%   Scalp data is simulated by projecting all components' signal activations 
%   through their respective, oriented source projections, and summing all projections together.

% Create [chan, samples, epochs] matrix
scalpdata = generate_scalpdata(c, lf, epochs);

%   Adds time stamps, channel names, etc.
epochs.marker = 'event 1';  % the epochs' time-locking event marker
epochs.prestim = 200;       % pre-stimulus period in ms. note that this
                            % only affects the time indicated in the final
                            % dataset; it is ignored during simulation.
                            % i.e., a simulated latency of 500 ms becomes a
                            % latency of 300 ms when a prestimulus period
                            % of 200 ms is later applied. you can use
                            % utl_shift_latency to shift all latencies in a
                            % class to fit the pre-stimulus period.

                            % Generate data in EEGlab format
EEG = utl_create_eeglabdataset(scalpdata, epochs, lf);

% Use eeglab to plot topos
pop_topoplot(EEG, 1, [100, 200, 250, 300, 350, 400, 500], '', [1 8]);

% Use eeglab to plot trials
pop_eegplot(EEG, 1, 1, 1);

%% Add ERP variability
%   add variability to the signal by indicating allowed random deviations or 
%   shifts as well as specific slopes. A deviation of 50 ms for our peak latency
%   allows this latency to vary +/- 50 ms between trials, following a normal
%   distribution with the indicated deviation being the six-sigma range.
erp.peakLatencyDv = 50; % applies to individual peaks
erp.peakAmplitudeDv = .2;
% erp.peakLatencyShift applies to all values equally

c.signal = {erp};
EEG = utl_create_eeglabdataset(generate_scalpdata(c, lf, epochs), ...
        epochs, lf);

% Change amplitude across epochs:
%   Indicating a slope results in a consistent change over time.
%   An amplitude of 1 and an amplitude slope of -.75, for example, results 
%   in the signal having a peak amplitude of 1 in the first epoch, and .25 
%   in the last.
erp.peakAmplitudeSlope = -.75;

c.signal = {erp};
EEG = utl_create_eeglabdataset(generate_scalpdata(c, lf, epochs), ...
        epochs, lf);

figure; pop_erpimage(EEG, 1, [25], [[]], 'Pz', 10, 1, {}, [], '', ...
        'yerplabel', '\muV', 'erp', 'on', 'cbar', 'on', ...
        'topo', {[25] EEG.chanlocs EEG.chaninfo});
pop_eegplot(EEG, 1, 1, 1);

% Plot extremes
plot_signal_fromclass(erp, epochs);

% Can also set peak probability which control whether it shows up on every trial

%% Add Noise
noise = struct( ...
        'type', 'noise', ...
        'color', 'brown', ...
        'amplitude', 1);
noise = utl_check_class(noise);

c.signal = {noise};
EEG = utl_create_eeglabdataset(generate_scalpdata(c, lf, epochs), ...
        epochs, lf);
pop_eegplot(EEG, 1, 1, 1);

% Combine ERP and noise
c.signal = {erp, noise};
EEG = utl_create_eeglabdataset(generate_scalpdata(c, lf, epochs), ...
        epochs, lf);
pop_eegplot(EEG, 1, 1, 1);

% Alternative using 3rd party uniform distribution instead of DSP toolbox Gaussian:
%   white-unif, brown-unif, etc. (white, pink, brown, blue, purple)

%% Simulate event-related spectral perturbations (ERSPs)
% Single frequency
ersp = struct( ...
        'type', 'ersp', ...
        'frequency', 20, ...
        'amplitude', .25, ...
        'modulation', 'none');
ersp = utl_check_class(ersp);
plot_signal_fromclass(ersp, epochs);

% Broadband (multiple frequencies)
%   maximum spectral power between 15 and 25 Hz, with transitions between 12-15
%   and 25-28 Hz. This is generated by filtering uniform white noise in the
%   indicated frequency band. In case a frequency band is indicated, phase
%   will be ignored.
ersp = struct( ...
        'type', 'ersp', ...
        'frequency', [12 15 25 28], ...
        'amplitude', .25, ...
        'modulation', 'none');

ersp = utl_check_class(ersp);

plot_signal_fromclass(ersp, epochs);

%% Add modulations:
% 1) single frequency burst with given peak (center) latency, width, taper
ersp.modulation = 'burst';
ersp.modLatency = 500;      % centre of the burst, in ms
ersp.modWidth = 100;        % width (half duration) of the burst, in ms
ersp.modTaper = 0.5;        % taper of the burst

ersp = utl_check_class(ersp);

plot_signal_fromclass(ersp, epochs);

% 2) inverse burst for event-related desynchronization
ersp.modulation = 'invburst';
ersp.modWidth = 300;
ersp.modMinRelAmplitude = 0.05;

plot_signal_fromclass(ersp, epochs);

% 3) Amplitude modulation according to phase of another
%   20 Hz base freq modulated with 2 Hz sine wave
ersp.frequency = 20;
ersp.modulation = 'ampmod';
ersp.modFrequency = 2;
ersp.modPhase = .25;

ersp = utl_check_class(ersp);

plot_signal_fromclass(ersp, epochs);
% Can use modMinRelAmplitude and phase options too

% Remove oscillation during baseline epoch
ersp.frequency = [12 15 25 28];
ersp.modPrestimPeriod = epochs.prestim;
ersp.modPrestimTaper = .5;

plot_signal_fromclass(ersp, epochs);

% Combine
c.signal = {erp, ersp};%, noise};
EEG = utl_create_eeglabdataset(generate_scalpdata(c, lf, epochs), ...
        epochs, lf);
pop_eegplot(EEG, 1, 1, 1);

%% Options using pre-generated data
% creating matrix of specific random data
s = rng(0, 'v5uniform');
randomdata = randn(epochs.n, epochs.srate*epochs.length/1000);
rng(s);

% adding generated random data to data class
data = struct();
data.data = randomdata;
data.index = {'e', ':'};
data.amplitude = .5;
data.amplitudeType = 'relative';

data = utl_check_class(data, 'type', 'data');

plot_signal_fromclass(data, epochs);

%% Auto-regressive models
% obtaining two data classes containing ARM-generated signals: two time
% series series with one unidirectional interaction between them, and a
% model order of 10
arm = arm_get_class_interacting(2, 10, 1, epochs, 1);

% getting two components with sources that are 100 mm apart
armsourcelocs = lf_get_source_spaced(lf, 2, 100);
armcomps = utl_create_component(armsourcelocs, arm, lf);

plot_component(armcomps, epochs, lf);

%% Random class generation for many components
% obtain 10 random source locations, at least 5 cm apart
sourcelocs  = lf_get_source_spaced(leadfield, 10, 50);

% obtain 64 random ERP activation classes. each class will have
% either 1, 2, or 3 peaks, each centred between 200 and 1000 ms,
% widths between 25 and 200 ms, and amplitudes between -1 and 1.
erp = erp_get_class_random([1:3], [200:1000], [25:200], ...
		[-1:.1:-.5, .5:.1:1], 'numClasses', 10);

% combining into brain components
c = utl_create_component(sourcelocs, erp, lf);

% plotting each component's projection and activation
plot_component(c, epochs, lf);

%% Source identification
% obtaining 64 random source locations, at least 2.5 cm apart
sourcelocs  = lf_get_source_spaced(leadfield, 64, 25);

% each source will get the same type of activation: brown coloured noise
signal      = struct('type', 'noise', 'color', 'brown', 'amplitude', 1);

% packing source locations and activations together into components
components = utl_create_component(sourcelocs, signal, lf);

% obtaining mixing and unmixing matrices
[w, winv]   = utl_get_icaweights(components, leadfield);