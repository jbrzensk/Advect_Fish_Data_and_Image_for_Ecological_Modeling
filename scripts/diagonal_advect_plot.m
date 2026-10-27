%% diagonal_advect_plot.m
% Compares final fish populations from three diagonal advection techniques:
%   1. Semi-Lagrangian Advection only              (test2_Diag_A)
%   2. Random Semi-Lagrangian Advection            (test2_Diag_RA)
%   3. Semi-Lagrangian Advection + Laplacian Diff  (test2_Diag_AD)
% AND plots them with their initial conditions and the true solution.
%
% Food gradient is diagonal (x+y), fish swim at 45 degrees.
% True solution advances with speed fish_speed/sqrt(2) in each direction.
%
% Loops over four CFL values (0.5, 0.9, 1.0, 1.25), loading the
% corresponding .mat files. Run the three test2_Diag scripts first.
%
% Brzenski, Jared
% 03/26/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clearvars;

addpath('../src')
addpath('../external/cmocean')

svgfilename = 'test2_plot1_fulla.svg';
filename    = 'test2_plot1_fulla.png';

% Two Column figure size (pixels), same as test1_plot1.png
imagesizevector = [0 0 1642 991];

width_a = 0.25;   % panel (a) width
width   = 0.22;
height  = 0.31;
%            1     2     3     4     5    cb
xlocs = [  0.075 0.39  0.68  0.39  0.68  0.91];
ylocs = [  0.33  0.51  0.51  0.09  0.09  0.09];

% Arrow of motion in panel (a), in data units [km]
xarrow = [540 615];
yarrow = [40 115];

%% Settings
cfl_targets = [0.5, 0.9, 1.0, 1.25];

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');
set(groot, ...
    'defaultAxesFontSize', 16, ...
    'defaultTextFontSize', 16, ...
    'defaultAxesFontName', 'Arial', ...
    'defaultTextFontName', 'Arial', ...
    'defaultLineLineWidth', 2.0);

labels = {'Advection (A)', 'Random Advect (RA)', 'Advect + Diffusion (AD)'};
%            1      2      3      4      5
sublabels = {'(a)', '(b)', '(c)', '(d)', '(e)'};

contour_colors  = {'b',           'r',          [0.1 0.7 0.1]};
contour_lstyles = {'-',           '--',         '-.'          };
contour_lwidths = [1.5,           1.5,          1.5           ];

titlefontsize     = 14;
minititlefontsize = 14;
labelfontsize     = 16;
cbfontsize        = 16;
legendfontsize    = 12;
sublabelfontsize  = 18;

%% =========================================================
%  Build diagonal initial conditions
% ==========================================================
days          = 2.5;
fish_speed_ms = 1.0;
fish_init     = 1.0;
food_init     = 1000;
day           = 24*60*60;
epsilon       = 10e-15;

x0_fish = 500;
y0_fish = 0;

m  = 200;   n  = 200;
a  = 0;     b  = 1000;
c  = -500.0; d = 500.0;
dx = (b-a)/m;
dy = (d-c)/n;
[X_ic, Y_ic] = meshgrid(a:dx:b, c:dy:d);

Fish_base        = fish_init * ones(n+1, m+1);
Fish_base        = zeroCornerMatrix(Fish_base);
Fish_base        = makeFood('humpfish', Fish_base, X_ic, Y_ic, x0_fish, y0_fish, 1, 0);
Fish_base(Fish_base < epsilon) = 0;

Food_base = ((Y_ic + 250) + X_ic) / 1000;
Food_base = Food_base .* food_init;

% Analytic solution at t = 2.5 days (45-degree diagonal)
distance        = fish_speed_ms * day * days / 1000 / sqrt(2);   % km
delta_x_final   = x0_fish + distance;
delta_y_final   = y0_fish + distance;
True_Fish_final = makeFood('humpfish', Fish_base, X_ic, Y_ic, delta_x_final, delta_y_final, 1, 0);
True_Fish_final(True_Fish_final < epsilon) = 0;

%% =========================================================
%  Create figure
% ==========================================================
fig1 = figure('Name', 'Diagonal Gradient', ...
              'Units', 'pixels', 'Position', imagesizevector);
colormap(fig1, flipud(gray));

global_clim_max = 0;
ax_all = gobjects(1, numel(cfl_targets));

%% =========================================================
%  Panel 1: Initial conditions (food + initial + analytic fish)
% ==========================================================
clim_fish   = max([Fish_base(:); True_Fish_final(:)]);
lev_min     = clim_fish * 0.03;
levels_fish = linspace(lev_min, clim_fish, 25);

ax_ic = axes('Position', [xlocs(1) ylocs(1) width_a height]);

contourf(ax_ic, X_ic, Y_ic, Food_base, 15, 'LineColor', 'none');
colormap(ax_ic, cmocean('haline'));
cb_ic = colorbar(ax_ic);
cb_ic.Location          = 'northoutside';
cb_ic.Label.String      = 'Habitat Suitability Field';
cb_ic.Label.Interpreter = 'latex';
cb_ic.TickLabelInterpreter = 'latex';
cb_ic.FontSize          = cbfontsize;
hold(ax_ic, 'on');

[~, hf] = contourf(ax_ic, X_ic, Y_ic, True_Fish_final, levels_fish, ...
                   'LineColor', 'k', 'LineWidth', 2);
[~, hi] = contourf(ax_ic, X_ic, Y_ic, Fish_base, levels_fish, ...
                   'LineColor', 'r', 'LineWidth', 2);

hold(ax_ic, 'off');

xlabel(ax_ic, 'x [km]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
ylabel(ax_ic, 'y [km]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
% Legend below the panel, as in horizontal_advect_plot.m
legend(ax_ic, [hf, hi], ...
       {sprintf('Analytic fish at $t=%.1f$d', days), 'Initial fish'}, ...
       'Interpreter', 'latex', 'Location', 'southoutside', 'FontSize', legendfontsize);
%xlim(ax_ic, [0 1000]);
%ylim(ax_ic, [-200 400]);
set(ax_ic, 'FontSize', labelfontsize, 'TickLabelInterpreter', 'latex');

% Equal x/y scale so the initial conditions are not distorted: set the axes
% height from its width and the data ranges, keeping the panel centred
axis(ax_ic, 'tight');
axpos  = get(ax_ic, 'Position');
figpos = get(fig1, 'Position');
h_eq   = axpos(3) * diff(ax_ic.YLim) / diff(ax_ic.XLim) * figpos(3) / figpos(4);
axpos  = [axpos(1), axpos(2) + (axpos(4) - h_eq)/2, axpos(3), h_eq];
set(ax_ic, 'Position', axpos);

% Diagonal arrow showing direction of motion, converted from data units to figure units
xfig  = axpos(1) + (xarrow - ax_ic.XLim(1)) / diff(ax_ic.XLim) * axpos(3);
yfig  = axpos(2) + (yarrow - ax_ic.YLim(1)) / diff(ax_ic.YLim) * axpos(4);
annotation('arrow', xfig, yfig);

% Panel label, above the colorbar label (normalized axes units)
text(ax_ic, -0.20, 1.29, ['\textbf{', sublabels{1}, '}'], 'Units', 'normalized', ...
     'Interpreter', 'latex', 'FontSize', sublabelfontsize, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

%% =========================================================
%  CFL loop — panels 2–5
% ==========================================================
for ci = 1:numel(cfl_targets)

    target_cfl = cfl_targets(ci);
    cfl_tag    = sprintf('cfl%03d', round(target_cfl * 100));
    cfl_str    = sprintf('CFL = %.2f', target_cfl);

    %% Load saved results for this CFL
    A  = load(['test2_Diag_A_',  cfl_tag, '.mat']);
    RA = load(['test2_Diag_RA_', cfl_tag, '.mat']);
    AD = load(['test2_Diag_AD_', cfl_tag, '.mat']);

    data  = {A, RA, AD};
    X     = A.X;
    Y     = A.Y;

    %% Error metrics
    fprintf('\n========== %s — Error Analysis ==========\n', cfl_str);
    fprintf('%-26s %10s %12s %10s %10s\n', 'Method', 'MSE', 'Correlation', 'E\_diss', 'E\_disp');
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

    %% Build CFL panel
    fields_2d = {A.Fish, RA.Fish, AD.Fish};
    clim_max  = max([A.Fish(:); RA.Fish(:); AD.Fish(:); A.True_Fish(:)]);
    n_levels  = 15;
    lev_min   = clim_max * 0.03;
    levels    = linspace(lev_min, clim_max, n_levels);

    global_clim_max = max(global_clim_max, clim_max);

    lct  = [xlocs(ci+1), ylocs(ci+1), width, height];
    ax1  = axes('Position', lct);
    ax_all(ci) = ax1;

    [~, h_true] = contourf(ax1, X, Y, A.True_Fish, levels, 'LineColor', 'none');
    clim(ax1, [0, clim_max]);
    hold(ax1, 'on');

    h = gobjects(3, 1);
    for k = 1:3
        [~, h(k)] = contour(ax1, X, Y, fields_2d{k}, levels, ...
                            'LineColor', contour_colors{k}, ...
                            'LineStyle', contour_lstyles{k}, ...
                            'LineWidth', contour_lwidths(k));
    end
    hold(ax1, 'off');

    if ci == 2
        legend(ax1, [h_true; h], ['True Solution (filled)', labels], ...
               'Location', 'northeast', 'FontSize', legendfontsize, 'Interpreter', 'latex');
    end
    grid on;

    if ci > 2
        xlabel(ax1, 'x [km]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
    end
    if ci == 1 || ci == 3
        ylabel(ax1, 'y [km]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
    else
        set(ax1, 'TickLength', [0 0])
    end

    xlim(ax1, [500 800]);
    ylim(ax1, [0   300]);
    title(ax1, cfl_str, 'Interpreter', 'latex', 'FontSize', minititlefontsize);

    % Panel label, top-left above the axes (level with the title)
    text(ax1, 0, 1.01, ['\textbf{', sublabels{ci+1}, '}'], 'Units', 'normalized', ...
         'Interpreter', 'latex', 'FontSize', sublabelfontsize, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

end % CFL loop

%% Shared colorbar — unify clim across all CFL tiles, attach at right
for ci = 1:numel(cfl_targets)
    clim(ax_all(ci), [0, 1]);
end
cb_ax = axes(fig1, 'Position', [xlocs(end) ylocs(end) 0.001 ylocs(2)+height-ylocs(end)], 'Visible', 'off');
clim(cb_ax, [0, 1]);
colormap(cb_ax, flipud(gray));
cb_shared = colorbar(cb_ax, 'Location', 'eastoutside', ...
                     'TickLabelInterpreter', 'latex', ...
                     'FontSize', cbfontsize);
cb_shared.Label.String      = 'Fish Concentration';
cb_shared.Label.Interpreter = 'latex';
cb_shared.FontSize          = cbfontsize;

print(gcf, svgfilename, '-dsvg');
exportgraphics(gcf, filename, 'Resolution', 300);
