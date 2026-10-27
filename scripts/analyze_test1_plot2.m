%% analyze_test1.m
% Compares final fish populations from three advection techniques:
%   1. Semi-Lagrangian Advection only              (test1_Linear_A)
%   2. Random Semi-Lagrangian Advection            (test1_Linear_RA)
%   3. Semi-Lagrangian Advection + Laplacian Diff  (test1_Linear_A_D)
%
% Loops over three CFL values (0.75, 1.0, 1.25), loading the corresponding
% .mat files and generating the same three plots for each CFL.
%
% Run the three test scripts first to generate the .mat save files, then
% run this script from the scripts/ directory.
%
% Brzenski, Jared
% 03/09/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clearvars;

%% Settings
% Save figure
basename = 'test1_plot2b';   % change to whatever name you want

cfl_targets = [0.5, 0.9, 1.0, 1.25];

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');
set(groot, ...
    'defaultAxesFontSize', 16, ...
    'defaultTextFontSize', 16, ...
    'defaultAxesFontName', 'Arial', ...
    'defaultTextFontName', 'Arial');

% Storage for Y-averaged profiles across all CFL runs (filled inside loop)
profiles_all = cell(numel(cfl_targets), 4);   % (cfl, method+truth)
x_vec_store  = [];

labels = {'Advection (A)', 'Random Advect (RA)', 'Advect + Diffusion (AD)'};
%            1      2      3      4
sublabels = {'(a)', '(b)', '(c)', '(d)'};

contour_colors  = {'b',           'r',          [0.1 0.7 0.1]};
contour_lstyles = {'-',           '--',         '-.'          };
contour_lwidths = [1.5,           1.5,          1.5           ];

line_colors  = {'b',  'r',  [0.1 0.7 0.1],  'k'};
line_lstyles = {'-',  '--', '-.',            ':'};
line_lwidths = [1.5,  1.5,  1.5,            2.5];

titlefontsize = 18;
minititlefontsize = 16;
labelfontsize = 16;
cbfontsize = 16;
legendfontsize = 14;
sublabelfontsize = 18;

%% Storage for shared colorbar (filled during loop)
global_clim_max = 0;
ax_all = gobjects(1, numel(cfl_targets));

%% CFL loop — one full set of plots per CFL value
for ci = 1:numel(cfl_targets)

    target_cfl = cfl_targets(ci);
    cfl_tag    = sprintf('cfl%03d', round(target_cfl * 100));
    cfl_str    = sprintf('CFL = %.2f', target_cfl);

    %% Load saved results for this CFL
    A  = load(['test1_Linear_A_',   cfl_tag, '.mat']);
    RA = load(['test1_Linear_RA_',  cfl_tag, '.mat']);
    AD = load(['test1_Linear_A_D_', cfl_tag, '.mat']);

    data  = {A, RA, AD};
    X     = A.X;
    Y     = A.Y;
    x_vec = X(1,:);   % 1D x-axis [km]

    %% =========================================================
    %  Error Metrics vs. True Solution
    % ==========================================================
    fprintf('\n========== %s — Error Analysis ==========\n', cfl_str);
    fprintf('%-26s %10s %12s %10s %10s\n', 'Method', 'MSE', 'Correlation', 'E_diss', 'E_disp');
    fprintf('%s\n', repmat('-', 1, 72));

    metrics = struct('mse', {}, 'rho', {}, 'e_diss', {}, 'e_disp', {});
    for k = 1:3
        d     = data{k};
        F_vec = d.Fish(:);
        T_vec = d.True_Fish(:);
        N     = numel(F_vec);

        mse    = (1/N) * sum((F_vec - T_vec).^2);
        rho    = corr(F_vec, T_vec);
        stdF   = std(F_vec);   stdT = std(T_vec);
        meanF  = mean(F_vec);  meanT = mean(T_vec);
        e_diss = (stdF - stdT)^2 + (meanF - meanT)^2;
        e_disp = 2*(1 - rho)*stdF*stdT;

        metrics(k).mse    = mse;
        metrics(k).rho    = rho;
        metrics(k).e_diss = e_diss;
        metrics(k).e_disp = e_disp;

        fprintf('%-26s %10.4f %12.4f %10.4f %10.4f\n', ...
                labels{k}, mse, rho, e_diss, e_disp);
    end
    fprintf('%s\n', repmat('-', 1, 72));
    [~, best] = min([metrics.mse]);
    fprintf('Lowest MSE: %s\n\n', labels{best});

    %% =========================================================
    %  Figure 1: 2D final fish populations — all on one panel
    %  True Solution as filled background; methods as contour lines
    % ==========================================================
    fields_2d = {A.Fish, RA.Fish, AD.Fish};
    clim_max  = max([A.Fish(:); RA.Fish(:); AD.Fish(:); A.True_Fish(:)]);
    % Focus contour levels on the non-zero range so all levels are visible
    n_levels  = 15;
    lev_min   = clim_max * 0.03;   % start at 3% of peak to skip the zero floor
    levels    = linspace(lev_min, clim_max, n_levels);

    %% Store Y-averaged profiles for combined figure after loop
    fields_1d = {A.Fish, RA.Fish, AD.Fish, A.True_Fish};
    for k = 1:4
        profiles_all{ci, k} = mean(fields_1d{k}, 1);
    end
    x_vec_store = x_vec;

end % CFL loop


%% =========================================================
%  Combined Figure: Y-averaged profiles, one panel per CFL
% ==========================================================
plot_labels = [labels, {'True Solution'}];

fig_profiles = figure('Name', 'Y-Averaged Profiles — All CFL', ...
                      'units', 'inches', 'Position', [1 1 10 9]);
tl2  = tiledlayout(2, 2, 'TileSpacing', 'compact');

for ci = 1:numel(cfl_targets)
    %subplot(2, 2, ci);
    ax1 = nexttile(tl2);
    hold on;
    for k = 1:4
        plot(x_vec_store, profiles_all{ci, k}, ...
             'Color',     line_colors{k}, ...
             'LineStyle', line_lstyles{k}, ...
             'LineWidth', line_lwidths(k));
    end
    hold off;

    if ci > 2
        xlabel(ax1, 'x [km]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
    else
       set(gca, 'XTick') 
    end

    if ci == 2 || ci == 4
       set(gca, 'TickLength', [0 0])
    end

    title(sprintf('CFL = %.2f', cfl_targets(ci)), 'Interpreter', 'latex', 'FontSize', minititlefontsize);
    grid on;
    xlim([400 700]);
    ylim([0  0.1])
    if ci == 1
        % Legend in its own row below the panels so it never covers data
        lgd = legend(ax1, plot_labels, 'Interpreter', 'latex', ...
                     'FontSize', legendfontsize, 'NumColumns', 2);
        lgd.Layout.Tile = 'south';
    end

    % Panel label, top-left above the axes
    text(ax1, -0.02, 1.04, ['\textbf{', sublabels{ci}, '}'], 'Units', 'normalized', ...
         'Interpreter', 'latex', 'FontSize', sublabelfontsize, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end
% One y-label shared by both rows (too long to repeat on each panel)
ylabel(tl2, 'Fish Concentration  (y-averaged)', 'Interpreter', 'latex', ...
       'FontSize', labelfontsize);
%sgtitle('Final Fish Populations: Y-Averaged Profiles', 'FontSize', titlefontsize,...
%           'Interpreter', 'latex','FontWeight', 'bold');

exportgraphics(fig_profiles, [basename '.png'], 'Resolution', 300);
exportgraphics(fig_profiles, [basename '.eps'], 'ContentType', 'vector');
print(fig_profiles, [basename '.svg'], '-dsvg');