%% cfl_shift_plot.m
% Merges the three CFL slice figures (Paper/Fig/cfl_lt_one, cfl_eq_one,
% cfl_gt_one, previously saved by hand from plotFishCompare.m) into one
% plot: the fish concentration along y = y_slice after linear advection
% (test 1 setup, A only) for a CFL below, equal to and above one, together
% with the analytic solution.
%
% Every run takes a whole number of steps to the same final time, so all
% three share one analytic solution. The final time is 3 days rounded up
% to a whole number of grid cells of travel, so the CFL = 1 run lands
% exactly on the grid.
%
% Run this script from the scripts/ directory.
%
% Brzenski, Jared
% 10/23/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clearvars;

addpath('../src')

svgfilename = 'cfl_shift_plot.svg';
filename    = 'cfl_shift_plot.png';

% Single Column figure size (pixels): half the width of the two-column
% advection figures, so the fonts print at the same size
imagesizevector = [0 0 821 616];

%% Settings
cfl_targets = [0.72, 1.0, 1.08];   % below, equal to and above one
y_slice     = 0;                   % y of the plotted slice [km]

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');
set(groot, ...
    'defaultAxesFontSize', 16, ...
    'defaultTextFontSize', 16, ...
    'defaultAxesFontName', 'Arial', ...
    'defaultTextFontName', 'Arial', ...
    'defaultLineLineWidth', 2.0);

%            CFL < 1   CFL = 1   CFL > 1        Analytic
line_colors  = {'b',      'r',      [0.1 0.7 0.1], 'k'};
line_lstyles = {'-',      '--',     '-.',          ':'};
line_lwidths = [1.5,      1.5,      1.5,           2.5];
line_markers = {'o',      's',      '^',           'none'};

titlefontsize  = 14;
labelfontsize  = 16;
legendfontsize = 12;

%% Fixed parameters (same as test1_Linear_A.m)
days          = 3;      % sim time [days]
fish_speed_ms = 1.0;    % fish speed [m/s]
eating_factor = 0.00;   % fraction of food eaten per dt
fish_init     = 1.0;    % fish biomass [g/m^3]
food_init     = 1000;   % food biomass [g/m^3]
u_max         = 0.;     % current speed [m/s]
v_max         = 0.0;
day           = 24*60*60;   % seconds in a day
epsilon       = 10e-15;     % tiny value
%
%% Food scenario parameters
x0         = 800;
y0         = 0;
amp        = food_init;
other_conc = 0.1 * food_init;
x0_fish    = 300;   % initial fish hump center [km]
%
%% Grid
m  = 100;   n  = 50;
a  = 0;     b  = 1000;
c  = -250.0; d = 250.0;
dx   = (b-a)/m;           % km
dy   = (d-c)/n;           % km
dx_m = dx * 1000.0;       % m
dy_m = dy * 1000.0;       % m
[X, Y] = meshgrid(a:dx:b, c:dy:d);
%
%% Initial conditions (reset each CFL run from these)
Fish_base        = fish_init * ones(n+1, m+1);
Fish_base        = zeroCornerMatrix(Fish_base);
Fish_base        = makeFood('humpfish', Fish_base, X, Y, x0_fish, 0, 1, 0);
Fish_base(Fish_base < epsilon) = 0;
%
Food_base = food_init * ones(n+1, m+1);
Food_base = makeFood('gradient', Food_base, X, Y, x0, y0, amp, other_conc);
Food_base = fliplr(Food_base);
Food_base = Food_base .* food_init;
%
mask = ones(size(Food_base));
%
% Current in m/day
current = zeros(n+1, m+1, 2);
current(:,:,1) = mask .* u_max .* (1 - (2.*Y/(2*d)).^2);
current(:,:,2) = mask .* v_max;
current = current .* day;   % m/s -> m/day

%% Final time: 3 days of travel, rounded up to whole grid cells
n_cells   = ceil(days * day * fish_speed_ms / dx_m);   % cells travelled
t_final_s = n_cells * dx_m / fish_speed_ms;            % [s]

[~, iy] = min(abs(Y(:,1) - y_slice));   % grid row of the slice
x_vec   = X(1,:);

%% =========================================================
%  CFL loop: run each case to t_final and keep the slice
% ==========================================================
slices     = zeros(numel(cfl_targets), m+1);
cfl_actual = zeros(1, numel(cfl_targets));

for ci = 1:numel(cfl_targets)

    % Whole number of steps to t_final, as close as possible to the target CFL
    iters = round(t_final_s * fish_speed_ms / (cfl_targets(ci) * dx_m));
    dt_s  = t_final_s / iters;              % seconds
    dt    = dt_s / day;                     % days
    cfl_actual(ci) = fish_speed_ms * dt_s / dx_m;

    fish_speed          = fish_speed_ms * day;   % m/day
    food_consumption_dt = eating_factor * dt;

    fprintf('CFL = %.3f  |  dt = %.4f h  |  iters = %d\n', ...
            cfl_actual(ci), dt * 24, iters);

    Fish = Fish_base .* mask;
    Food = Food_base .* mask;
    for it = 1:iters
        Fish = Fish .* mask;
        [Fish, Food] = SemiLagrangianAdvectFish( Fish, ...
                                                 Food, ...
                                                 current, ...
                                                 food_consumption_dt, ...
                                                 dt, ...
                                                 dx_m, ...
                                                 dy_m, ...
                                                 fish_speed, ...
                                                 mask, ...
                                                 m, ...
                                                 n);
    end
    slices(ci,:) = Fish(iy,:);

end % CFL loop

%% Analytic solution at t_final, on a fine grid along the slice
% makeFood zeros the first and last rows, so evaluate on three rows and
% keep the middle one
x_fine     = linspace(a, b, 2001);
x_true     = x0_fish + t_final_s * fish_speed_ms / 1000;   % hump center [km]
[Xf, Yf]   = meshgrid(x_fine, Y(iy,1) + [-dy 0 dy]);
True_fine  = makeFood('humpfish', zeros(size(Xf)), Xf, Yf, x_true, 0, 1, 0);
True_slice = True_fine(2,:);
True_slice(True_slice < epsilon) = 0;

%% =========================================================
%  Figure: slices for every CFL plus the analytic solution
% ==========================================================
cfl_labels = arrayfun(@(c) sprintf('CFL = %.2f', c), cfl_actual, 'UniformOutput', false);

fig1 = figure('Name', 'CFL Shift Comparison', ...
              'Units', 'pixels', 'Position', imagesizevector);
ax = axes(fig1);
hold(ax, 'on');

% Analytic first so the numerical solutions are drawn on top
h_true = plot(ax, x_fine, True_slice, ...
              'Color',     line_colors{end}, ...
              'LineStyle', line_lstyles{end}, ...
              'LineWidth', line_lwidths(end));

h = gobjects(numel(cfl_targets), 1);
for ci = 1:numel(cfl_targets)
    h(ci) = plot(ax, x_vec, slices(ci,:), ...
                 'Color',     line_colors{ci}, ...
                 'LineStyle', line_lstyles{ci}, ...
                 'LineWidth', line_lwidths(ci), ...
                 'Marker',    line_markers{ci}, ...
                 'MarkerSize', 7);
end
hold(ax, 'off');

set(ax, 'FontSize', labelfontsize, 'TickLabelInterpreter', 'latex');
xlim(ax, [400 700]);
ylim(ax, [0 1.1]);
xlabel(ax, 'x [km]', 'Interpreter', 'latex', 'FontSize', labelfontsize);
ylabel(ax, 'Fish Concentration', 'Interpreter', 'latex', 'FontSize', labelfontsize);
%title(ax, sprintf('Slice at y = %g km, t = %.2f days', Y(iy,1), t_final_s/day), ...
%      'Interpreter', 'latex', 'FontSize', titlefontsize);

legend(ax, [h; h_true], [cfl_labels, {'Analytic'}], 'Interpreter', 'latex', ...
       'Location', 'northwest', 'FontSize', legendfontsize);
grid(ax, 'on');
box(ax, 'on');

print(fig1, svgfilename, '-dsvg');
exportgraphics(fig1, filename, 'Resolution', 300);
