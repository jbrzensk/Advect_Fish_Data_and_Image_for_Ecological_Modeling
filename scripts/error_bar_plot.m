%% error_bar_plot.m
% Rebuilds the combined error bar chart (previously assembled in Inkscape as
% Combined_Errors_bar_chart.svg) directly in MATLAB:
%   Row 1: Horizontal advection  (test1_Linear_* .mat files)
%   Row 2: Diagonal advection    (test2_Diag_*   .mat files)
%   Columns: E_diss, E_disp and MSE versus the true solution, with one bar
%   group per method and one bar per CFL value.
% All six panels share the same y-axis limits.
%
% Run the test1_Linear and test2_Diag scripts first to generate the .mat
% save files, then run this script from the scripts/ directory.
%
% Brzenski, Jared
% 10/14/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clearvars;

svgfilename = 'error_bar_plot.svg';
filename    = 'error_bar_plot.png';

imagesizevector = [1 1 15 10];   % [inches]

%% Settings
cfl_targets = [0.5, 0.9, 1.0, 1.25];

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');
set(groot, ...
    'defaultAxesFontSize', 16, ...
    'defaultTextFontSize', 16, ...
    'defaultAxesFontName', 'Arial', ...
    'defaultTextFontName', 'Arial');

labels = {'Advection (A)', 'Random Advect (RA)', 'Advect + Diffusion (AD)'};

% One row per test case: row title and .mat file prefixes for A, RA, A+D
case_titles = {'Horizontal Advection', 'Diagonal Advection'};
case_files  = {{'test1_Linear_A_', 'test1_Linear_RA_', 'test1_Linear_A_D_'}, ...
               {'test2_Diag_A_',   'test2_Diag_RA_',   'test2_Diag_AD_'   }};

metric_titles  = {'$\varepsilon$ dissipative', '$\varepsilon$ dispersive', 'MSE'};
metric_ylabels = {'$\mathbf{E}_{diss}$', '$\mathbf{E}_{disp}$', '\textbf{MSE}'};

%            1      2      3      4      5      6
sublabels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

ystep = 0.5e-3;   % shared y-limit is rounded up to a multiple of this

titlefontsize     = 18;
rowtitlefontsize  = 18;
minititlefontsize = 16;
labelfontsize     = 16;
legendfontsize    = 14;
sublabelfontsize  = 18;

%% =========================================================
%  Error metrics vs. true solution for every case, CFL and method
%  metric_data{tc, mp} is (n_cfl x 3 methods); mp = 1 E_diss, 2 E_disp, 3 MSE
% ==========================================================
metric_data = cell(numel(case_titles), 3);
for tc = 1:numel(case_titles)
    ediss_all = zeros(numel(cfl_targets), 3);
    edisp_all = zeros(numel(cfl_targets), 3);
    mse_all   = zeros(numel(cfl_targets), 3);

    for ci = 1:numel(cfl_targets)
        cfl_tag = sprintf('cfl%03d', round(cfl_targets(ci) * 100));
        for k = 1:3
            d     = load([case_files{tc}{k}, cfl_tag, '.mat'], 'Fish', 'True_Fish');
            F_vec = d.Fish(:);
            T_vec = d.True_Fish(:);
            N     = numel(F_vec);

            rho   = corr(F_vec, T_vec);
            stdF  = std(F_vec);   stdT = std(T_vec);
            meanF = mean(F_vec);  meanT = mean(T_vec);

            ediss_all(ci,k) = (stdF - stdT)^2 + (meanF - meanT)^2;
            edisp_all(ci,k) = 2*(1 - rho)*stdF*stdT;
            mse_all(ci,k)   = (1/N) * sum((F_vec - T_vec).^2);
        end
    end
    metric_data(tc,:) = {ediss_all, edisp_all, mse_all};
end

% Shared y-axis for all panels
ymax = max(cellfun(@(m) max(m(:)), metric_data(:)));
ymax = ceil(ymax / ystep) * ystep;

%% =========================================================
%  Figure: 2 rows (cases) x 3 columns (metrics)
%  Each panel: methods on x-axis, one bar group per CFL value
% ==========================================================
cfl_labels  = arrayfun(@(c) sprintf('CFL=%.2f', c), cfl_targets, 'UniformOutput', false);
labels_bold = cellfun(@(s) ['\textbf{' s '}'], labels, 'UniformOutput', false);

fig1 = figure('Name', 'Error Metrics — All CFL', ...
              'Units', 'inches', 'Position', imagesizevector);
% 'normal' padding leaves room for the rotated tick labels at the left edge
tl_outer = tiledlayout(fig1, numel(case_titles), 1, 'TileSpacing', 'compact', 'Padding', 'normal');
%title(tl_outer, 'Error by Type', 'Interpreter', 'latex', 'FontSize', titlefontsize);

for tc = 1:numel(case_titles)
    % Nested layout per row so each row gets its own title
    tl_row = tiledlayout(tl_outer, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    tl_row.Layout.Tile = tc;
    title(tl_row, case_titles{tc}, 'Interpreter', 'latex', 'FontSize', rowtitlefontsize);

    for mp = 1:3
        ax = nexttile(tl_row);
        bar(ax, metric_data{tc,mp}');          % size: (3 methods) x (4 CFL)
        set(ax, 'XTick', 1:3, 'XTickLabel', labels_bold, ...
                'TickLabelInterpreter', 'latex', 'FontSize', labelfontsize);
        xtickangle(ax, 20);
        if tc < numel(case_titles)
            set(ax, 'XTickLabel', []);         % method names on bottom row only
        end
        ylim(ax, [0 ymax]);
        ylabel(ax, metric_ylabels{mp}, 'Interpreter', 'latex', 'FontSize', labelfontsize);
        title(ax,  metric_titles{mp},  'Interpreter', 'latex', 'FontSize', minititlefontsize);
        grid(ax, 'on');
        if mp == 1
            legend(ax, cfl_labels, 'Interpreter', 'latex', ...
                   'Location', 'northwest', 'FontSize', legendfontsize);
        end

        % Panel label, top-left above the axes (left of the x10^-3 exponent)
        text(ax, -0.02, 1.04, ['\textbf{', sublabels{(tc-1)*3 + mp}, '}'], 'Units', 'normalized', ...
             'Interpreter', 'latex', 'FontSize', sublabelfontsize, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end
end

print(gcf, svgfilename, '-dsvg');
exportgraphics(gcf, filename, 'Resolution', 300);
