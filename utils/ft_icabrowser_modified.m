function [rej_comp] = ft_icabrowser(SBJ, cfg, comp)

% ICA component viewer and GUI
%
% loads in comp structure from FieldTrip ft_componentanalysis
% presents a GUI interface showin the power spectrum, variance over time
% and the topography of the components, as well as the possibility to save
% a PDF, view the timecourse and toggle components to be rejected vs kept.
% when done, will create a file with the components to be rejected
%
% CONFIGURATION NEEDED:
% cfg.path         where pdfs will be saves
% cfg.prefix       prefix of the pdf files
% cfg.layout       layout of the topo view
%
% OPTIONAL CONFIGURATION:
% cfg.colormap      colormap for topo
% cfg.inputfile
% cfg.outputfile    will contain indices of all components to reject
%
% original written by Thomas Pfeffer
% adapted by Jonathan Daume and Anne Urai
% University Medical Center Hamburg-Eppendorf, 2015

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(cfg, 'inputfile'), load(cfg.inputfile); end
assert(exist('comp', 'var') > 0, 'could not find a comp structure');

% setup
var_data = cat(2,comp.trial{:});
var_time = (1:(size(var_data,2)))/comp.fsample;
% only do the fft on a subset of trials, saves time
fft_data = cat(2,comp.trial{1:5:end});
% preallocate rejected components
rej_comp = zeros(size(comp.label,1),1);

subpl = 1;
l = 0;
cnt = 1;
il = 0;
logar = 1; % default logarithmic or linear scale for powspctrm

% set path
path = cfg.path;
if ~exist(path, 'dir'),
    mkdir path;
end
prefix = cfg.prefix;

% to save time redoing this for each topo
cfglay = keepfields(cfg, {'layout'});
comp2 = comp;
comp2.label = comp.topolabel;
lay = ft_prepare_layout(cfglay, comp2);

cfgtopo = [];
cfgtopo.layout    = lay;     % specify the layout file that should be used for plotting
cfgtopo.comment   = 'no';
cfgtopo.highlight = 'off';
cfgtopo.marker    = 'off';
cfgtopo.style     = 'straight';
if isfield(cfg, 'colormap'),  cfgtopo.colormap  = cfg.colormap;  end

err = 0;
manpos = [0.1 0.1 0.8 0.8]; % figure position, can be updated later

% ------------------------------------------------
% COMPUTE LATENCY FOR 2s-WINDOWS
% ------------------------------------------------

slen = floor(2*comp.fsample);
smax = floor(size(var_data,2)/slen);
comp_var  = nan(subpl, smax); % preallocate
comp_time = nan(1, smax);     % preallocate
for s = 1 : smax
  comp_time(s) = mean(var_time(1,(s-1)*slen+1:s*slen));
end
for x = 1: numel(comp.label)       
        il = il + 1;
        i = (cnt-1)*subpl+il;
        
            % keep manual screen position - better in dual monitor settings
fig_name = [SBJ comp.label(il,1)]
fig_name = strcat(fig_name(1), fig_name(2), '_fig1');
fig_name = char(fig_name);
f = figure('units','normalized','outerposition', manpos, 'Name', fig_name, 'Visible', 'off');
l = l + 1;

        
        % ------------------------------------------------
        % COMPUTE VARIANCE FOR 2s-WINDOWS
        % ------------------------------------------------
        
        for s = 1 : smax
            comp_var(i,s)=var(var_data(i,(s-1)*slen+1:s*slen));
        end
        
        % ------------------------------------------------
        % COMPUTE POWER SPECTRUM
        % ------------------------------------------------
        smo = 50;
        steps = 10;
        Fs = comp.fsample;
        N = floor(size(fft_data,2));
        xdft = fft(fft_data(i,:));
        xdft = xdft(1:N/2+1);
        psdx = (1/(Fs*N)).*abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        
        j = 1;
        k = 1;
        while j < length(psdx)-smo
            smoothed(k)=mean(psdx(j:j+smo));
            j = j + steps;
            k = k + 1;
        end
        
        freq = linspace(0,Fs/2,size(smoothed,2));
        strt = find(freq > 2,1,'first');
        stp  = find(freq < 200,1,'last');
        
        % ------------------------------------------------
        % PLOT POWER SPECTRUM
        % ------------------------------------------------
        subcomp{1}{il} = subplot('Position', [0.1, 0.1, 0.4, 0.4]);
        if logar
            plot(freq(strt:stp),log10(smoothed(strt:stp)));
            ylabel('(dB/Hz)');
        else
            plot(freq(strt:stp),smoothed(strt:stp));
            ylabel('T^2/Hz');
        end
        set(gca,'TickDir','out','XTick',0:25:200)
        xlabel('Frequency (Hz)'); grid on;
        axis tight;
        
        % ------------------------------------------------
        % PLOT VARIANCE OVER TIME
        % ------------------------------------------------
        subcomp{2}{il} = subplot('Position', [0.1, 0.6, 0.4, 0.4]);
        scatter(comp_time,comp_var(i,:),'k.');
        xlabel('Time (s)'); ylabel('Variance');
        axis tight; set(gca, 'tickdir', 'out');
        
        % ------------------------------------------------
        % PLOT COMPONENT TOPOGRAPHY
        % ------------------------------------------------
        subcomp{3}{il} = subplot('Position', [0.55, 0.1, 0.4, 0.9]);
        cfgtopo.component = i;       % specify the component(s) that should be plotted
        ft_topoplotIC(cfgtopo, comp);
        
        if mod(i,subpl)==0 || i == 80
            
            pos = [0.76 0.73 0.075 0.035; ...
                0.76 0.51 0.075 0.035; ...
                0.76 0.29 0.075 0.035; ...
                0.76 0.07 0.075 0.035];
        end
        fig_path = [cfg.path fig_name '.png'];
        saveas(f, fig_path);
    end
end


  