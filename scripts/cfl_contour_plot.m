%% cfl_contour_plot.m
% Rebuilds the CFL sweep error figure (previously assembled in Inkscape as
% cfl_tests_all.svg) directly in MATLAB:
%   Row 1: Horizontal advection  (test3: CFL_A/RA/AD_error.mat)
%   Row 2: Diagonal advection    (test5: CFL_Diag_A/RA/AD_error.mat)
%   Columns: one method each. Filled contours of total error versus maximum
%   speed and time step, with CFL contours (black) and the boundary of the
%   lowest 25% error region (red) overlaid.
% All six panels share the same color scale.
%
% Run the test3 and test5 scripts first to generate the .mat save files,
% then run this script from the scripts/ directory.
%
% Brzenski, Jared
% 10/14/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clearvars;

svgfilename = 'cfl_contour_plota.svg';
filename    = 'cfl_contour_plota.png';

imagesizevector = [1 1 16 11];   % [inches]

%% Settings
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');
set(groot, ...
    'defaultAxesFontSize', 16, ...
    'defaultTextFontSize', 16, ...
    'defaultAxesFontName', 'Arial', ...
    'defaultTextFontName', 'Arial');

% One row per test case: row title and .mat files for A, RA, A+D
case_titles = {'Horizontal Advection', 'Diagonal Advection'};
case_files  = {{'CFL_A_error',      'CFL_RA_error',      'CFL_AD_error'     }, ...
               {'CFL_Diag_A_error', 'CFL_Diag_RA_error', 'CFL_Diag_AD_error'}};
meth_labels = {'Advection (A)', 'Random Advection (RA)', 'Advection $+$ Diffusion (AD)'};

%            1      2      3      4      5      6
sublabels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

num_contours = 10;
cfl_levels   = linspace(0.1, 2, 20);
clim_max     = 4;     % E_tot above this saturates (horizontal RA peaks near 5.5)

rowtitlefontsize  = 18;
minititlefontsize = 16;
labelfontsize     = 16;
clabelfontsize    = 10;
sublabelfontsize  = 18;

%% Load data
data = cell(numel(case_titles), 3);
for tc = 1:numel(case_titles)
    for k = 1:3
        data{tc,k} = load([case_files{tc}{k}, '.mat']);
    end
end

% Shared color limits across all panels
clim_min    = min(cellfun(@(d) min(d.scaled_E_total(:)), data(:)));
shared_clim = [clim_min, clim_max];

%% =========================================================
%  Figure: 2 rows (cases) x 3 columns (methods)
%  Each panel: total error (filled) + CFL overlay
% ==========================================================
fig1 = figure('Name', 'Total Error with CFL Contours — All Methods', ...
              'Units', 'inches', 'Position', imagesizevector);
tl_outer = tiledlayout(fig1, numel(case_titles), 1, 'TileSpacing', 'compact', 'Padding', 'normal');

for tc = 1:numel(case_titles)
    % Nested layout per row so each row gets its own title
    tl_row = tiledlayout(tl_outer, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    tl_row.Layout.Tile = tc;
    title(tl_row, case_titles{tc}, 'Interpreter', 'latex', 'FontSize', rowtitlefontsize);

    for k = 1:3
        d  = data{tc,k};
        Xp = d.X_plt;
        Yp = d.Y_plt;

        ax = nexttile(tl_row);
        contourf(ax, Xp, Yp, d.scaled_E_total', num_contours, 'LineStyle', 'none');
        colormap(ax, parula);
        clim(ax, shared_clim);
        hold(ax, 'on');

        % 25% lowest-error region boundary
        sorted_vals = sort(d.scaled_E_total(:));
        thresh      = sorted_vals(round(0.25 * numel(sorted_vals)));
        contour(ax, Xp, Yp, d.scaled_E_total', [thresh thresh], ...
                'r', 'LineWidth', 2.5, 'DisplayName', 'min 25\% region');

        % CFL contours
        [C, h] = contour(ax, Xp, Yp, d.cfl', cfl_levels, ...
                         'k', 'LineWidth', 1.0, 'DisplayName', 'CFL');
        clabel(C, h, 'FontSize', clabelfontsize, 'LabelSpacing', 400);
        hold(ax, 'off');

        pbaspect(ax, [1 1 1]);
        if tc == numel(case_titles)
            xlabel(ax, '$S_{\max}$ [m/s]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
        end
        if k == 1
            ylabel(ax, '$\Delta t$ [hr]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
        end
        title(ax, meth_labels{k}, 'Interpreter', 'latex', 'FontSize', minititlefontsize);
        set(ax, 'FontSize', labelfontsize, 'TickLabelInterpreter', 'latex');

        % Panel label, top-left above the axes
        text(ax, -0.02, 1.04, ['\textbf{', sublabels{(tc-1)*3 + k}, '}'], 'Units', 'normalized', ...
             'Interpreter', 'latex', 'FontSize', sublabelfontsize, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end

    % One colorbar per row; every panel uses the same color limits
    cb = colorbar(ax);
    cb.Layout.Tile = 'east';
    ylabel(cb, '$\epsilon_{tot}\ (\times 10^{-3})$', ...
           'Interpreter', 'latex', 'FontSize', labelfontsize);
    set(cb, 'TickLabelInterpreter', 'latex', 'FontSize', labelfontsize);
end

print(gcf, svgfilename, '-dsvg');
exportgraphics(gcf, filename, 'Resolution', 300);
