%% Test 4: Kappa sensitivity for Advection + Diffusion
%  Fish follow food — linear advection with varying diffusion coefficient
%
% Runs 5 cases at fixed CFL = 1.0:
%   kappa = 0      → advection only (no laplacian2D), matches test1_Linear_A
%   kappa = 250    → advection + weak diffusion
%   kappa = 500    → advection + moderate diffusion
%   kappa = 750    → advection + strong diffusion
%   kappa = 1000   → advection + max diffusion (matches test1_Linear_A_D)
%
% All results compared against the true (exact) solution.
% Saves: test4_kappa_AD.mat
%
% BY: Jared Brzenski
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clearvars

addpath('../src')

savefilename = 'test4_kappa_AD.mat';

%% Sweep parameters
kappa_tries = [0, 250, 500, 750, 1000];   % [m^2/s]; 0 = advection only
n_cases     = numel(kappa_tries);

%% Fixed simulation parameters
target_cfl    = 1.0;
days          = 3;
fish_speed_ms = 1.0;       % [m/s]
eating_factor = 0.00;
fish_init     = 1.0;
food_init     = 1000;
u_max         = 0.0;
v_max         = 0.0;
day           = 24*60*60;  % seconds per day
epsilon       = 10e-15;

%% Food / fish scenario
x0         = 800;
y0         = 0;
amp        = food_init;
other_conc = 0.1 * food_init;
x0_fish    = 300;          % initial fish hump center [km]

%% Grid
m    = 100;   n    = 50;
a    = 0;     b    = 1000;
c    = -250.0; d   = 250.0;
dx   = (b-a)/m;
dy   = (d-c)/n;
dx_m = dx * 1000.0;        % [m]
dy_m = dy * 1000.0;        % [m]

[X, Y] = meshgrid(a:dx:b, c:dy:d);

vol   = dx * dy * 0.1;           % cell volume [km^3], 100 m thick
vol_m = vol * (1000^2 * 100);    % cell volume [m^3]

%% Time discretization from CFL
fish_speed = fish_speed_ms * day;                       % [m/day]
dt_s       = target_cfl * dx_m / fish_speed_ms;         % [s]
dt         = dt_s / day;                                 % [days]
num_hours  = dt * 24;
iters      = ceil(days * 24 / num_hours);
food_consumption_dt = eating_factor * dt;

fprintf('CFL = %.2f  |  dt = %.4f days (%.4f h)  |  iters = %d\n', ...
        target_cfl, dt, num_hours, iters);

%% Base initial conditions (shared; reset per kappa)
Fish_base = fish_init * ones(n+1, m+1);
Fish_base = zeroCornerMatrix(Fish_base);
Fish_base = makeFood('humpfish', Fish_base, X, Y, x0_fish, 0, 1, 0);
Fish_base(Fish_base < epsilon) = 0;

Food_base = food_init * ones(n+1, m+1);
Food_base = makeFood('gradient', Food_base, X, Y, x0, y0, amp, other_conc);
Food_base = fliplr(Food_base);
Food_base = Food_base .* food_init;

mask = ones(size(Food_base));

current_base = zeros(n+1, m+1, 2);
current_base(:,:,1) = mask .* u_max .* (1 - (2.*Y/(2*d)).^2);
current_base(:,:,2) = mask .* v_max;
current_base = current_base .* day;   % m/s -> m/day

%% Storage
Fish_all  = zeros(n+1, m+1, n_cases);   % final fish field per kappa
E_total   = zeros(1, n_cases);
E_diss    = zeros(1, n_cases);
E_disp    = zeros(1, n_cases);
% Fish mass history: mass_total(ki,k) [g] at time mass_time(k) [days]
mass_time  = (0:iters) * dt;
mass_total = zeros(n_cases, iters+1);

%% =========================================================
%  Kappa loop
%% =========================================================
for ki = 1:n_cases
    kappa = kappa_tries(ki);
    kappa_dt = kappa * dt_s;    % [m^2]; 0 when kappa=0

    if kappa == 0
        fprintf('\n--- kappa = 0  (advection only) ---\n');
    else
        fprintf('\n--- kappa = %d m^2/s  |  kappa_dt = %.2f m^2 ---\n', kappa, kappa_dt);
    end

    %% Reset fields
    Fish    = Fish_base;
    Food    = Food_base;
    current = current_base;
    current(:,:,1) = squeeze(current(:,:,1)) .* mask;
    current(:,:,2) = squeeze(current(:,:,2)) .* mask;
    Food = Food .* mask;
    Fish = Fish .* mask;

    %% Initial true solution
    delta_x   = x0_fish + (0 * dt_s * fish_speed_ms / 1000);
    True_Fish = makeFood('humpfish', Fish_base, X, Y, delta_x, 0, 1, 0);
    True_Fish(True_Fish < epsilon) = 0;

    mass_total(ki,1) = sum(Fish, 'all') * vol_m;

    %% Iteration loop
    for it = 1:iters
        Fish = Fish .* mask;
        Fish_init = Fish;

        [Fish, Food, ~] = SemiLagrangianAdvectFish( Fish, ...
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

        % Apply diffusion only when kappa > 0
        if kappa > 0
            Fish = laplacian2D(Fish, dx_m, dy_m, kappa_dt);
        end

        % Advance true solution
        delta_x   = x0_fish + (it * (fish_speed * dt) / 1000);
        True_Fish = makeFood('humpfish', Fish_init, X, Y, delta_x, 0, 1, 0);
        True_Fish(True_Fish < epsilon) = 0;

        % Record total fish mass after this step
        mass_total(ki,it+1) = sum(Fish, 'all') * vol_m;
    end

    %% Error calculations
    F_vec = Fish(:);
    A_vec = True_Fish(:);

    E_total(ki) = 1/(m*n) * sum((F_vec - A_vec).^2);
    corrFA      = corr(F_vec, A_vec);
    stdF        = std(F_vec);   stdA  = std(A_vec);
    meanF       = mean(F_vec);  meanA = mean(A_vec);
    E_diss(ki)  = (stdF - stdA)^2 + (meanF - meanA)^2;
    E_disp(ki)  = 2*(1 - corrFA)*stdF*stdA;

    fprintf('MSE=%.4e  E_diss=%.4e  E_disp=%.4e\n', E_total(ki), E_diss(ki), E_disp(ki));
    fprintf('Fish mass relative change: %+.3e\n', ...
            (mass_total(ki,end) - mass_total(ki,1)) / mass_total(ki,1));

    Fish_all(:,:,ki) = Fish;
end

%% Save
save(savefilename, 'Fish_all', 'True_Fish', 'X', 'Y', ...
     'kappa_tries', 'E_total', 'E_diss', 'E_disp', ...
     'dx_m', 'dy_m', 'dt', 'dt_s', 'target_cfl', ...
     'mass_time', 'mass_total');
fprintf('\nSaved: %s\n', savefilename);

%% =========================================================
%  Plotting — publication quality
%% =========================================================
labelfont    = 12;
titlefont    = 13;
lineweight   = 1.5;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');

% Build display labels
case_labels = cell(1, n_cases);
for ki = 1:n_cases
    if kappa_tries(ki) == 0
        case_labels{ki} = 'Advection Only';
    else
        case_labels{ki} = sprintf('$\\kappa = %d$ m$^2$/s', kappa_tries(ki));
    end
end

% Consistent color cycle
cmap = lines(n_cases);

%% Shared color limits across all panels
clim_max = max(True_Fish(:)) * 1.05;
clim_min = 0;

%% Figure 1: 2D fish distributions — True + 5 kappa cases
fig1 = figure('Units', 'inches', 'Position', [1 1 13 5]);
t1   = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'tight');

% Panel 1: True Solution
ax0 = nexttile;
contourf(X, Y, True_Fish, 10, 'LineStyle', 'none');
colormap(ax0, flipud(gray));
clim([clim_min clim_max]);
colorbar(ax0);
axis tight; pbaspect([2 1 1]);
xlabel('$x$ [km]', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$y$ [km]', 'Interpreter', 'latex', 'FontSize', labelfont);
title('True Solution', 'Interpreter', 'latex', 'FontSize', titlefont);
set(ax0, 'FontSize', labelfont, 'TickLabelInterpreter', 'latex');

% Panels 2–6: kappa cases
panel_letters = {'(b)','(c)','(d)','(e)','(f)'};
for ki = 1:n_cases
    ax = nexttile;
    contourf(X, Y, Fish_all(:,:,ki), 10, 'LineStyle', 'none');
    colormap(ax, flipud(gray));
    clim([clim_min clim_max]);
    colorbar(ax);
    axis tight; pbaspect([2 1 1]);
    xlabel('$x$ [km]', 'Interpreter', 'latex', 'FontSize', labelfont);
    ylabel('$y$ [km]', 'Interpreter', 'latex', 'FontSize', labelfont);
    if kappa_tries(ki) == 0
        ttl = sprintf('%s Advection Only', panel_letters{ki});
    else
        ttl = sprintf('%s $\\kappa = %d$ m$^2$/s', panel_letters{ki}, kappa_tries(ki));
    end
    title(ttl, 'Interpreter', 'latex', 'FontSize', titlefont);
    set(ax, 'FontSize', labelfont, 'TickLabelInterpreter', 'latex');
end

title(t1, 'Fish Density at $t = 3$ days  (CFL $= 1.0$)', ...
      'Interpreter', 'latex', 'FontSize', titlefont+1);

%% Figure 2: Y-averaged profiles
fig2 = figure('Units', 'inches', 'Position', [1 1 7 5]);

% True solution
y_avg_true = mean(True_Fish, 1);
x_vec      = a:dx:b;
plot(x_vec, y_avg_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True Solution');
hold on;

for ki = 1:n_cases
    y_avg = mean(Fish_all(:,:,ki), 1);
    plot(x_vec, y_avg, 'LineWidth', lineweight, ...
         'Color', cmap(ki,:), 'DisplayName', case_labels{ki});
end
hold off;

axis tight;
xlabel('$x$ [km]', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\bar{C}$ [g m$^{-3}$]', 'Interpreter', 'latex', 'FontSize', labelfont);
title('$y$-Averaged Fish Density at $t = 3$ days', ...
      'Interpreter', 'latex', 'FontSize', titlefont);
legend('Location', 'northwest', 'FontSize', labelfont-1);
set(gca, 'FontSize', labelfont, 'TickLabelInterpreter', 'latex');
grid on; box on;

%% Figure 3: Error metrics vs kappa
fig3 = figure('Units', 'inches', 'Position', [1 1 7 5]);

bar_data = [E_diss; E_disp; E_total]';   % [n_cases x 3]
b = bar(kappa_tries, bar_data, 'grouped');
b(1).DisplayName = '$E_{diss}$';
b(2).DisplayName = '$E_{disp}$';
b(3).DisplayName = '$\epsilon_{tot}$';

xlabel('$\kappa$ [m$^2$ s$^{-1}$]', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('Error',                       'Interpreter', 'latex', 'FontSize', labelfont);
title('Error vs.\ Diffusion Coefficient (CFL $= 1.0$)', ...
      'Interpreter', 'latex', 'FontSize', titlefont);
legend('Interpreter', 'latex', 'Location', 'northwest', 'FontSize', labelfont-1);
set(gca, 'FontSize', labelfont, 'TickLabelInterpreter', 'latex');
xticks(kappa_tries);
grid on; box on;
