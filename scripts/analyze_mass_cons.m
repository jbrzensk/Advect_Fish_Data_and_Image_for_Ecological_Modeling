%% analyze_mass_cons.m
% Fish mass conservation across all test cases.
%
% Loads mass_time / mass_total from the .mat files written by:
%   test1_Linear_{A,A_D,RA}.m     per-CFL time series (linear)
%   test2_Diag_{A,AD,RA}.m        per-CFL time series (diagonal)
%   test4_kappa_AD.m              per-kappa time series
%   test3_CFL_{A,AD,RA}_error.m   speed x dt sweeps (linear)
%   test5_diag_CFL_{A,AD,RA}_error.m  speed x dt sweeps (diagonal)
%
% The drift is (M(t) - M0)/M0. A conservative scheme gives values at the
% round-off level, ~1e-16, so every plot shows the drift in multiples of
% machine epsilon (eps = 2.2e-16). A value of 3 means 3*eps, i.e. the
% scheme conserves mass to the last few bits of a double.
%
% Figures:
%   1. Drift vs time for test1, test2, test4 (shared y-scale)
%   2. Worst drift over time for every run in the test3/test5 sweeps
%   3. Summary: worst drift of every run, grouped by experiment
%
% Files that are missing, or that predate mass tracking, are skipped with
% a warning. Run from the scripts/ directory.
%
% Brzenski, Jared
% 10/16/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clearvars;

%% Settings
cfl_targets = [0.5, 0.9, 1.0, 1.25];
save_figs   = 1;          % 1 = write PNGs of each figure

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');

titlefontsize = 13;
labelfontsize = 12;
tickfontsize  = 11;

% Chart colors
c_surface = [252 252 251] / 255;
c_ink     = [ 11  11  11] / 255;
c_ink2    = [ 82  81  78] / 255;   % secondary text, reference lines
c_muted   = [137 135 129] / 255;   % axes
c_grid    = [225 224 217] / 255;
c_series  = [ 42 120 214;          % blue    (CFL 0.50 / kappa 0)
             235 104  52;          % orange  (CFL 0.90 / kappa 250)
              27 175 122;          % aqua    (CFL 1.00 / kappa 500)
             237 161   0;          % yellow  (CFL 1.25 / kappa 750)
             232 123 164] / 255;   % magenta (kappa 1000)
% Sequential blue ramp (light = near zero, dark = large)
blue_ramp = [205 226 251; 158 197 244; 109 167 236; 57 135 229; ...
              42 120 214;  37 106 191;  24  79 149; 13  54 107] / 255;
cmap_seq  = interp1(linspace(0, 1, size(blue_ramp,1)), blue_ramp, linspace(0, 1, 256));

ep = eps;   % 2.2204e-16

%% Experiment definitions
% Single-run experiments: one file per CFL
single_prefix = {'test1_Linear_A_', 'test1_Linear_A_D_', 'test1_Linear_RA_', ...
                 'test2_Diag_A_',   'test2_Diag_AD_',    'test2_Diag_RA_'};
single_label  = {'Linear: Advection (A)', 'Linear: A + Diffusion', 'Linear: Random Advect (RA)', ...
                 'Diagonal: Advection (A)', 'Diagonal: A + Diffusion', 'Diagonal: Random Advect (RA)'};

% Sweep experiments: one file holding a speed x dt grid of runs
sweep_file  = {'CFL_A_error.mat',      'CFL_AD_error.mat',      'CFL_RA_error.mat', ...
               'CFL_Diag_A_error.mat', 'CFL_Diag_AD_error.mat', 'CFL_Diag_RA_error.mat'};
sweep_label = {'Linear sweep: A', 'Linear sweep: A + D', 'Linear sweep: RA', ...
               'Diagonal sweep: A', 'Diagonal sweep: A + D', 'Diagonal sweep: RA'};

% Summary collection: worst |drift|/eps of each run, per experiment
sum_label = {};
sum_vals  = {};

%% Load single-run experiments
n_single = numel(single_prefix);
n_cfl    = numel(cfl_targets);
ts_t     = cell(n_single, n_cfl);   % time [days]
ts_d     = cell(n_single, n_cfl);   % drift [eps]

for g = 1:n_single
    vals = [];
    for ci = 1:n_cfl
        cfl_tag = sprintf('cfl%03d', round(cfl_targets(ci) * 100));
        S = load_mass([single_prefix{g}, cfl_tag, '.mat']);
        if isempty(S), continue; end
        ts_t{g,ci} = S.mass_time;
        ts_d{g,ci} = (S.mass_total - S.mass_total(1)) / S.mass_total(1) / ep;
        vals(end+1) = max(abs(ts_d{g,ci}));  %#ok<SAGROW>
    end
    if ~isempty(vals)
        sum_label{end+1} = single_label{g};  %#ok<SAGROW>
        sum_vals{end+1}  = vals;             %#ok<SAGROW>
    end
end

%% Load test4 (kappa sweep, shared time vector)
K = load_mass('test4_kappa_AD.mat', {'kappa_tries'});
if ~isempty(K)
    k_d = (K.mass_total - K.mass_total(:,1)) ./ K.mass_total(:,1) / ep;
    sum_label{end+1} = 'Linear: kappa sweep (CFL 1)';
    sum_vals{end+1}  = max(abs(k_d), [], 2)';
end

%% Load sweep experiments
n_sweep  = numel(sweep_file);
sw       = cell(1, n_sweep);
for g = 1:n_sweep
    S = load_mass(sweep_file{g}, {'X_plt', 'Y_plt', 'cfl'});
    if isempty(S), continue; end
    % Worst |drift| over the whole run, for every (speed, dt) pair
    S.max_d = cellfun(@(M) max(abs(M - M(1))) / M(1), S.mass_total) / ep;
    sw{g} = S;
    sum_label{end+1} = sweep_label{g};  %#ok<SAGROW>
    sum_vals{end+1}  = S.max_d(:)';     %#ok<SAGROW>
end

if isempty(sum_vals)
    error(['No .mat files with mass_time/mass_total found. ' ...
           'Rerun the test scripts to generate them.']);
end

%% Console summary
fprintf('\nFish mass conservation, drift = (M(t) - M0)/M0, eps = %.3e\n', ep);
fprintf('%-32s %6s %14s %14s %14s\n', 'Experiment', 'Runs', 'max |drift|', ...
        'max [eps]', 'median [eps]');
for k = 1:numel(sum_vals)
    fprintf('%-32s %6d %14.3e %14.2f %14.2f\n', sum_label{k}, numel(sum_vals{k}), ...
            max(sum_vals{k}) * ep, max(sum_vals{k}), median(sum_vals{k}));
end
worst_all = max(cellfun(@max, sum_vals));
fprintf('Worst drift over all runs: %.3e  (%.2f eps)\n\n', worst_all * ep, worst_all);

%% =========================================================
%  Figure 1: drift vs time, test1 / test2 / test4
%% =========================================================
all_d = [ts_d(:); {[]}];
if ~isempty(K), all_d{end+1} = k_d(:)'; end
d_max = max([cellfun(@(d) max([abs(d), 0]), all_d); 0]);
y_lim = max(4, ceil(d_max + 0.5)) * [-1 1];   % shared, symmetric, at least +/-4 eps

fig1 = figure('Name', 'Mass drift vs time', 'Color', c_surface, ...
              'Units', 'inches', 'Position', [1 1 12 9.5]);
tl1  = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, 'Fish mass drift $(M(t)-M_0)/M_0$ in multiples of machine $\varepsilon$', ...
      'Interpreter', 'latex', 'FontSize', titlefontsize + 1);

h_cfl = gobjects(1, n_cfl);
for g = 1:n_single
    ax = nexttile(tl1, g);
    style_axes(ax, c_surface, c_muted, c_grid, tickfontsize);
    yline(ax, 0, 'Color', c_muted, 'LineWidth', 1);
    for ci = 1:n_cfl
        if isempty(ts_t{g,ci}), continue; end
        h = plot(ax, ts_t{g,ci}, ts_d{g,ci}, '-', 'Color', c_series(ci,:), 'LineWidth', 1.5);
        h_cfl(ci) = h;
    end
    ylim(ax, y_lim);
    xlim(ax, [0 3.1]);
    title(ax, single_label{g}, 'Interpreter', 'latex', 'FontSize', titlefontsize);
    if mod(g-1, 3) == 0
        ylabel(ax, 'drift [$\varepsilon$]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
    end
    xlabel(ax, 'time [days]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
end
ok = isgraphics(h_cfl);
if any(ok)
    cfl_names = arrayfun(@(c) sprintf('CFL = %.2f', c), cfl_targets, 'UniformOutput', false);
    lgd = legend(h_cfl(ok), cfl_names(ok), 'Orientation', 'horizontal', ...
                 'FontSize', labelfontsize, 'Box', 'off', 'TextColor', c_ink);
    lgd.Layout.Tile = 'north';
end

% test4 spans the bottom row
ax = nexttile(tl1, 7, [1 3]);
style_axes(ax, c_surface, c_muted, c_grid, tickfontsize);
yline(ax, 0, 'Color', c_muted, 'LineWidth', 1);
if ~isempty(K)
    h_k = gobjects(1, numel(K.kappa_tries));
    for ki = 1:numel(K.kappa_tries)
        h_k(ki) = plot(ax, K.mass_time, k_d(ki,:), '-', ...
                       'Color', c_series(ki,:), 'LineWidth', 1.5);
    end
    k_names = arrayfun(@(k) sprintf('$\\kappa = %d$', k), K.kappa_tries, 'UniformOutput', false);
    legend(ax, h_k, k_names, 'Location', 'eastoutside', 'FontSize', labelfontsize, ...
           'Box', 'off', 'TextColor', c_ink);
else
    text(ax, 0.5, 0.5, 'test4\_kappa\_AD.mat not available', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'Color', c_ink2);
end
ylim(ax, y_lim);
xlim(ax, [0 3.1]);
title(ax, 'Linear: A + Diffusion, $\kappa$ sweep at CFL = 1 [m$^2$/s]', ...
      'Interpreter', 'latex', 'FontSize', titlefontsize);
ylabel(ax, 'drift [$\varepsilon$]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
xlabel(ax, 'time [days]', 'Interpreter', 'latex', 'FontSize', labelfontsize);

%% =========================================================
%  Figure 2: worst drift over time for every sweep run
%% =========================================================
have_sw = ~cellfun(@isempty, sw);
if any(have_sw)
    c_max = max(cellfun(@(S) max(S.max_d(:)), sw(have_sw)));
    c_max = max(1, ceil(c_max));

    fig2 = figure('Name', 'Mass drift, CFL sweeps', 'Color', c_surface, ...
                  'Units', 'inches', 'Position', [1 1 12 8]);
    tl2  = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl2, 'Worst fish mass drift over each run, $\max_t |M(t)-M_0|/M_0$ [$\varepsilon$]', ...
          'Interpreter', 'latex', 'FontSize', titlefontsize + 1);
    colormap(fig2, cmap_seq);

    for g = 1:n_sweep
        ax = nexttile(tl2, g);
        if ~have_sw(g)
            axis(ax, 'off');
            text(ax, 0.5, 0.5, [strrep(sweep_file{g}, '_', '\_') ' not available'], ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center', 'Color', c_ink2);
            continue;
        end
        S = sw{g};
        speed = S.X_plt(1,:);
        dt_hr = S.Y_plt(:,1);
        imagesc(ax, speed, dt_hr, S.max_d');
        axis(ax, 'xy');
        hold(ax, 'on');
        % CFL = 1 for orientation
        if numel(speed) > 1 && numel(dt_hr) > 1
            [Cc, hc] = contour(ax, S.X_plt, S.Y_plt, S.cfl', [1 1], ...
                               'Color', c_ink2, 'LineWidth', 1);
            clabel(Cc, hc, 'FontSize', tickfontsize - 1, 'Color', c_ink2, ...
                   'LabelSpacing', 400);
        end
        hold(ax, 'off');
        clim(ax, [0 c_max]);
        set(ax, 'FontSize', tickfontsize, 'TickLabelInterpreter', 'latex', ...
                'XColor', c_muted, 'YColor', c_muted, 'Layer', 'top', 'Box', 'off');
        pbaspect(ax, [1 1 1]);
        title(ax, sweep_label{g}, 'Interpreter', 'latex', 'FontSize', titlefontsize);
        xlabel(ax, '$S_{\max}$ [m/s]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
        ylabel(ax, '$\Delta t$ [hr]',  'Interpreter', 'latex', 'FontSize', labelfontsize);
    end
    cb = colorbar;
    cb.Layout.Tile = 'east';
    set(cb, 'TickLabelInterpreter', 'latex', 'FontSize', tickfontsize, 'Color', c_muted);
    ylabel(cb, 'drift [$\varepsilon$]', 'Interpreter', 'latex', 'FontSize', labelfontsize, ...
           'Color', c_ink);
end

%% =========================================================
%  Figure 3: summary, every run grouped by experiment
%% =========================================================
n_exp = numel(sum_vals);
fig3  = figure('Name', 'Mass drift summary', 'Color', c_surface, ...
               'Units', 'inches', 'Position', [1 1 9 0.45*n_exp + 1.8]);
ax = axes(fig3);
style_axes(ax, c_surface, c_muted, c_grid, tickfontsize);
ax.YGrid = 'off';

rng(1);   % fixed jitter, same picture every run
x_hi = max(4, ceil(worst_all + 0.5));
for k = 1:n_exp
    v  = sum_vals{k};
    yk = n_exp - k + 1;    % first experiment on top
    jit = (rand(size(v)) - 0.5) * 0.5;
    scatter(ax, v, yk + jit, 22, c_series(1,:), 'filled', ...
            'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none');
    % Worst run in this experiment
    plot(ax, max(v) * [1 1], yk + [-0.32 0.32], '-', 'Color', c_ink, 'LineWidth', 2);
    text(ax, max(v) + 0.012 * x_hi, yk, sprintf('%.1f', max(v)), ...
         'Color', c_ink2, 'FontSize', tickfontsize - 1, 'VerticalAlignment', 'middle', ...
         'Interpreter', 'latex');
end
xlim(ax, [0 x_hi * 1.06]);
ylim(ax, [0.4 n_exp + 0.6]);
set(ax, 'YTick', 1:n_exp, 'YTickLabel', flip(sum_label), 'TickLength', [0 0], ...
        'YColor', c_ink2);
xlabel(ax, 'worst drift over run, $\max_t |M(t)-M_0|/M_0$ [$\varepsilon$]', ...
       'Interpreter', 'latex', 'FontSize', labelfontsize);
% Summing ~1e4 grid values in double precision can drift by up to ~1e4 eps,
% so anything below 1e-12 relative is round-off, not a leak in the scheme.
if worst_all * ep < 1e-12
    ttl = 'Fish mass is conserved to round-off';
else
    ttl = 'Fish mass drifts beyond round-off';
end
title(ax, sprintf('%s: worst run is %.1f$\\varepsilon$ (%.1e)', ttl, worst_all, worst_all * ep), ...
      'Interpreter', 'latex', 'FontSize', titlefontsize);
subtitle(ax, 'one dot per run; black bar = worst run in each experiment', ...
         'Interpreter', 'latex', 'FontSize', tickfontsize, 'Color', c_ink2);

%% Save figures
if save_figs
    exportgraphics(fig1, 'mass_cons_timeseries.png', 'Resolution', 200);
    if exist('fig2', 'var'), exportgraphics(fig2, 'mass_cons_sweeps.png', 'Resolution', 200); end
    exportgraphics(fig3, 'mass_cons_summary.png', 'Resolution', 200);
    fprintf('Saved mass_cons_*.png\n');
end

%% =========================================================
%  Local functions
%% =========================================================
function S = load_mass(fname, extra_vars)
% Load mass_time / mass_total (plus extra_vars) from fname.
% Returns [] with a warning if the file or the mass variables are missing.
    if nargin < 2, extra_vars = {}; end
    S = [];
    if ~isfile(fname)
        warning('analyze_mass_cons:missing', 'Skipping %s (file not found).', fname);
        return;
    end
    vars = who('-file', fname);
    if ~all(ismember({'mass_time', 'mass_total'}, vars))
        warning('analyze_mass_cons:nomass', ...
                'Skipping %s (no mass data, rerun its test script).', fname);
        return;
    end
    S = load(fname, 'mass_time', 'mass_total', extra_vars{:});
end

function style_axes(ax, c_surface, c_muted, c_grid, fs)
% Recessive axes: hairline solid grid, muted axis color, no box.
    hold(ax, 'on');
    set(ax, 'Color', c_surface, 'XColor', c_muted, 'YColor', c_muted, ...
            'GridColor', c_grid, 'GridAlpha', 1, 'GridLineStyle', '-', ...
            'XGrid', 'on', 'YGrid', 'on', 'Box', 'off', ...
            'FontSize', fs, 'TickLabelInterpreter', 'latex');
end
