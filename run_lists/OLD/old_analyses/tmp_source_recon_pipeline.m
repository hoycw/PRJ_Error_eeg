

ft_dir = [app_dir 'fieldtrip/'];

% Load Head Model (vol)
load([ft_dir 'template/headmodel/standard_bem.mat'])
headmodel = vol;

% Load Elec postitions
elec = ft_read_sens([ft_dir 'template/electrode/standard_1020.elc']);
% !!! will need to select my 64 electrodes
% !!! check if I need to realign these!
elec_align = elec;

% Load T1 MRI
mri = ft_read_mri([ft_dir 'template/anatomy/single_subj_T1.nii']);

% Load Source Model (sourcemodel)
load([ft_dir 'template/sourcemodel/standard_sourcemodel3d4mm.mat']);
% figure
% plot3(sourcemodel.pos(:,1), sourcemodel.pos(:,2), sourcemodel.pos(:,3), '.')
% axis equal
% axis vis3d
% grid on
% ft_plot_mesh(sourcemodel);


%% Create Forward Solution: Leadfield
cfg            = [];
cfg.headmodel  = headmodel;
cfg.elec       = elec_align;
cfg.grid       = sourcemodel;
%cfg.singleshell.batchsize = 5000; % speeds up the computation
leadfield      = ft_prepare_leadfield(cfg);
% cfg.channel = 'all'; (or subset)

% plot the leadfield for a few representative locations: points around
% z-axis with increasing z values

plotpos   = [];
positions = [];
n = size(leadfield.pos,1);
p = 1;
for i = 1:n
  if leadfield.pos(i,1)==-0.2 && leadfield.pos(i,2)==tmp(28)
      plotpos(p) = i;
      positions(p,:) = leadfield.pos(i,:);
      p=p+1;
  end
end
figure;
for i=1:20
    subplot(4,5,i);
    if ~isempty(leadfield.leadfield{plotpos(i)})
        ft_plot_topo3d(leadfield.cfg.elec.chanpos,leadfield.leadfield{plotpos(i)}(:,3));
        %view([0 0]);
    end
end

%% Load and Preprocess EEG Data


%% Inverse Solution:
cfg               = [];
cfg.method        = 'mne';
cfg.grid          = leadfield;
cfg.headmodel     = headmodel;
cfg.mne.prewhiten = 'yes';
cfg.mne.lambda    = 3;
cfg.mne.scalesourcecov = 'yes';
sourceFC          = ft_sourceanalysis(cfg, tlckFC);
sourceFIC         = ft_sourceanalysis(cfg, tlckFIC);

