%% horizontal_advect_plot.m
% Compares final fish populations from three advection techniques:
%   1. Semi-Lagrangian Advection only              (test1_Linear_A)
%   2. Random Semi-Lagrangian Advection            (test1_Linear_RA)
%   3. Semi-Lagrangian Advection + Laplacian Diff  (test1_Linear_A_D)
% AND plots them with thier initial conditions and the true solution.
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

addpath('../src')
addpath('../external/cmocean')

svgfilename = 'test1_plot1_fulla.svg';
filename = 'test1_plot1_fulla.png';
% Single Column
imagesizevector = [1 1 4 4];
% Two Column figure size (pixels), same as test1_plot1.png
imagesizevector = [0 0 1642 991];

width_a = 0.25;   % panel (a) is narrower so its colorbar label clears panel (b)
width = 0.22;
height = 0.31;
%           1     2    3    4    5   cb
xlocs = [  0.075 0.39 0.68 0.39 0.68 0.91];
ylocs = [  0.33 0.51 0.51 0.09 0.09 0.09];

% Arrow of motion in panel (a), in data units [km]
xarrow = [360 490];
yarrow = [0 0];

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

% Storage for Y-averaged profiles across all CFL runs (filled inside loop)
profiles_all = cell(numel(cfl_targets), 4);   % (cfl, method+truth)
x_vec_store  = [];

labels = {'Advection (A)', 'Random Advect (RA)', 'Advect + Diffusion (AD)'};
%            1      2      3      4      5
sublabels = {'(a)', '(b)', '(c)', '(d)', '(e)'};

contour_colors  = {'b',           'r',          [0.1 0.7 0.1]};
contour_lstyles = {'-',           '--',         '-.'          };
contour_lwidths = [1.5,           1.5,          1.5           ];

line_colors  = {'b',  'r',  [0.1 0.7 0.1],  'k'};
line_lstyles = {'-',  '--', '-.',            ':'};
line_lwidths = [1.5,  1.5,  1.5,            2.5];

titlefontsize = 14;
minititlefontsize = 14;
labelfontsize = 16;
cbfontsize = 16;
legendfontsize = 12;
sublabelfontsize = 18;

%% Create combined 2D figure before loop (panels filled inside)
fig1 = figure('Name', 'Horizontal Gradient', ...
              'Units', 'pixels', 'Position', imagesizevector);
%tl1  = tiledlayout(2, 3, 'TileSpacing', 'compact');
colormap(fig1, flipud(gray));
%title(fig1, 'Fish Concentration [$g/m^{3}$], methods versus analytic', ...
%      'Interpreter', 'latex', 'FontSize', titlefontsize);

%% Storage for shared colorbar (filled during loop)
global_clim_max = 0;
ax_all = gobjects(1, numel(cfl_targets));

%% Plot initial conditions
[Food_base, Fish_base, True_Fish_final, X, Y] = horizontal_initial_conditions();
% Tiled image #1
clim_fish   = max([Fish_base(:); True_Fish_final(:)]);
lev_min     = clim_fish * 0.03;
levels_fish = linspace(lev_min, clim_fish, 25);

ax = axes('Position',[xlocs(1) ylocs(1) width_a height]);

% Food gradient magnitude as background (gray colormap)
contourf(ax, X, Y, Food_base, 15, 'LineColor', 'none');
colormap(ax, cmocean('haline'));   % keep background colormap
cb = colorbar(ax);
cb.Location = 'northoutside';
cb.Label.String      = 'Habitat Suitability Field';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = cbfontsize;
hold(ax, 'on');

% Analytic fish at t=3d — filled contour (blue)
[~, hf] = contourf(ax, X, Y, True_Fish_final, levels_fish, ...
                  'LineColor', 'k', 'LineWidth', 2);


% Initial fish — contour lines only (red dashed)
[~, hi] = contourf(ax, X, Y, Fish_base, levels_fish, ...
                  'LineColor', 'r', 'LineWidth', 2);

hold(ax, 'off');

xlabel(ax, 'x [km]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
ylabel(ax, 'y [km]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
%title('Initial Fish, Analytic Solution & Food Concentration','Interpreter', 'latex', 'FontSize', fs_title);
% Legend below the panel; inside the short panel it would cover the fish
legend(ax, [hf, hi], ...
       {'Analytic fish at $t=3$d', 'Initial fish'}, ...
       'Interpreter', 'latex', 'Location', 'southoutside', 'FontSize', legendfontsize);
set(ax, 'FontSize', labelfontsize, 'TickLabelInterpreter', 'latex');

% Equal x/y scale so the initial conditions are not distorted: set the axes
% height from its width and the data ranges, keeping the panel centred
axis(ax, 'tight');
axpos  = get(ax, 'Position');
figpos = get(fig1, 'Position');
h_eq   = axpos(3) * diff(ax.YLim) / diff(ax.XLim) * figpos(3) / figpos(4);
axpos  = [axpos(1), axpos(2) + (axpos(4) - h_eq)/2, axpos(3), h_eq];
set(ax, 'Position', axpos);

% Arrow of motion, converted from data units to figure units
xfig  = axpos(1) + (xarrow - ax.XLim(1)) / diff(ax.XLim) * axpos(3);
yfig  = axpos(2) + (yarrow - ax.YLim(1)) / diff(ax.YLim) * axpos(4);
annotation('arrow', xfig, yfig);

% Panel label, above the colorbar label (normalized axes units)
text(ax, -0.20, 1.58, ['\textbf{', sublabels{1}, '}'], 'Units', 'normalized', ...
     'Interpreter', 'latex', 'FontSize', sublabelfontsize, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

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

    % Tiles 3-6

    fields_2d = {A.Fish, RA.Fish, AD.Fish};
    clim_max  = max([A.Fish(:); RA.Fish(:); AD.Fish(:); A.True_Fish(:)]);
    % Focus contour levels on the non-zero range so all levels are visible
    n_levels  = 15;
    lev_min   = clim_max * 0.03;   % start at 3% of peak to skip the zero floor
    levels    = linspace(lev_min, clim_max, n_levels);

    global_clim_max = max(global_clim_max, clim_max);
    % ax1 = nexttile(tl1);
    %ax1 = subplot(2,3,locations(ci) );
    lct = [xlocs(ci+1), ylocs(ci+1), width, height];
    ax1 = axes('Position', lct);
    ax_all(ci) = ax1;
    [~, h_true] = contourf(ax1, X, Y, A.True_Fish, levels, 'LineColor', 'none');
    clim(ax1, [0, clim_max]);
    hold(ax1, 'on');

    h = gobjects(3,1);
    % All three methods as contour lines
    for k = 1:3
        [~, h(k)] = contour(ax1, X, Y, fields_2d{k}, levels, ...
                            'LineColor',  contour_colors{k}, ...
                            'LineStyle',  contour_lstyles{k}, ...
                            'LineWidth',  contour_lwidths(k));
    end
    hold(ax1, 'off');
    
    if ci == 1
        legend(ax1, [h_true; h], ['True Solution (filled)', labels], ...
               'Location', 'northeast', 'FontSize', legendfontsize, 'Interpreter','latex');
    end
    grid on;

    if ci > 2
        xlabel(ax1, 'x [km]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
    else
       set(gca, 'XTick') 
    end

    if ci == 1 || ci == 3
        ylabel(ax1, 'y [km]', 'Interpreter', 'latex', 'FontSize',labelfontsize);
    else
       set(gca, 'TickLength', [0 0])
    end

    xlim(ax1, [360 660]);
    ylim(ax1, [-150 150]);
    title(ax1, cfl_str, 'Interpreter','latex', 'Fontsize', minititlefontsize);

    % Panel label, top-left above the axes (level with the title)
    text(ax1, 0, 1.01, ['\textbf{', sublabels{ci+1}, '}'], 'Units', 'normalized', ...
         'Interpreter', 'latex', 'FontSize', sublabelfontsize, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

    %% Store Y-averaged profiles for combined figure after loop
    fields_1d = {A.Fish, RA.Fish, AD.Fish, A.True_Fish};
    for k = 1:4
        profiles_all{ci, k} = mean(fields_1d{k}, 1);
    end
    x_vec_store = x_vec;

end % CFL loop

%% Shared colorbar — unify clim across all tiles, add one bar at the bottom
for ci = 1:numel(cfl_targets)
    clim(ax_all(ci), [0, 1]);
end
% Attach colorbar to a hidden axes so it spans the full figure width
cb_ax = axes(fig1, 'Position', [xlocs(end) ylocs(end) 0.001 ylocs(2)+height-ylocs(end)], 'Visible', 'off');
clim(cb_ax, [0, 1]);
colormap(cb_ax, flipud(gray));
cb_shared = colorbar(cb_ax, 'Location', 'eastoutside', ...
                     'TickLabelInterpreter', 'latex', ...
                     'FontSize', cbfontsize);
cb_shared.Label.String         = 'Fish Concentration';
cb_shared.Label.Interpreter    = 'latex';
cb_shared.FontSize             = cbfontsize;

print(gcf, svgfilename, '-dsvg');
exportgraphics(gcf, filename, 'Resolution', 300);
