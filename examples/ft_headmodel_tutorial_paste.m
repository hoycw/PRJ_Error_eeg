%% Create Head Model
% Load MRI
mri = ft_read_mri(mri_fname);

% Re-slice for homogeneous voxel size
%ft_volumereslice if necessary

% Align MRI to head coordinate system
%ft_determine_coordsys
%ft_volumerealign

% Segmentation
% FreeSurfer

% Prepare mesh
cfg = [];
cfg.tissue      = {'brain','skull','scalp'};
cfg.numvertices = [3000 2000 1000];
bnd = ft_prepare_mesh(cfg,mri_seg);

% Head model
cfg = [];
cfg.method = 'dipoli'; % or 'openmeeg', 'bemcp', etc
vol = ft_prepare_headmodel(cfg, bnd);

% Plot
figure;
ft_plot_mesh(vol.bnd(3),'facecolor','none'); % Scalp
figure;
ft_plot_mesh(vol.bnd(2),'facecolor','none'); % Skull
figure;
ft_plot_mesh(vol.bnd(1),'facecolor','none'); % Brain
% Combined
figure;
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(vol.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4]);

% Align electrodes
elec = ft_read_sens('standard_1020.elc');

nas=mri.hdr.fiducial.mri.nas;
lpa=mri.hdr.fiducial.mri.lpa;
rpa=mri.hdr.fiducial.mri.rpa;

transm=mri.transform;

nas=ft_warp_apply(transm,nas, 'homogenous');
lpa=ft_warp_apply(transm,lpa, 'homogenous');
rpa=ft_warp_apply(transm,rpa, 'homogenous');

% create a structure similar to a template set of electrodes
fid.elecpos       = [nas; lpa; rpa];       % ctf-coordinates of fiducials
fid.label         = {'Nz','LPA','RPA'};    % same labels as in elec
fid.unit          = 'mm';                  % same units as mri

% alignment
cfg               = [];
cfg.method        = 'fiducial';
cfg.target        = fid;                   % see above
cfg.elec          = elec;
cfg.fiducial      = {'Nz', 'LPA', 'RPA'};  % labels of fiducials in fid and in elec
elec_aligned      = ft_electroderealign(cfg);

figure;
ft_plot_sens(elec_aligned,'style','sk');
hold on;
ft_plot_mesh(vol.bnd(1),'facealpha', 0.85, 'edgecolor', 'none', 'facecolor', [0.65 0.65 0.65]); %scalp


% Save
% save(vol_fname,'-v7.3','vol');

%% Source model
