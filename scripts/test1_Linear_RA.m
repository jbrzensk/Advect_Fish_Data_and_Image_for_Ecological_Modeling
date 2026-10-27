%% Test_1, linear Random Advection of fish — CFL sweep
%  Fish follow food
%
% Test Case 1: Linear advection for 3 days.
% No consumption, no speed adjustment, just
% if Conc_n > Conc_m, swim towards Cell_n at max_speed
%
% Runs three simulations at CFL = 0.75, 1.0, 1.25 and saves a separate
% .mat file for each.  dt is derived from the target CFL each run.
%
% BY: Jared Brzenski
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clearvars

addpath('../src')

%% CFL targets — dt is computed from each
cfl_targets = [0.5, 0.9, 1.0, 1.25];

%% Fixed parameters (same for every CFL run)
days          = 3;      % sim time [days]
fish_speed_ms = 1.0;    % fish speed [m/s]  — converted to m/day inside loop
eating_factor = 0.00;   % fraction of food eaten per dt
fish_init     = 1.0;    % fish biomass [g/m^3]
food_init     = 1000;   % food biomass [g/m^3]
u_max         = 0.;     % current speed [m/s]
v_max         = 0.0;
day           = 24*60*60;   % seconds in a day
epsilon       = 10e-15;     % tiny value
%
make_plot  = 0;   % set to 1 to enable per-step plotting (slower)
make_video = 0;
do_smooth  = 0;
%
%% Food scenario parameters
x0         = 800;
y0         = 0;
amp        = food_init;
other_conc = 0.1 * food_init;
x0_fish    = 300;   % initial fish hump center [km]
%
%% Grid (fixed for all runs)
m  = 100;   n  = 50;
a  = 0;     b  = 1000;
c  = -250.0; d = 250.0;
dx   = (b-a)/m;           % km
dy   = (d-c)/n;           % km
dx_m = dx * 1000.0;       % m
dy_m = dy * 1000.0;       % m
%
vol   = dx * dy * 0.1;
vol_m = vol * (1000^2 * 100);
length_u = (m+1)*n;
length_v = m*(n+1);
length_c = (m+2)*(n+2);
[X, Y] = meshgrid(a:dx:b, c:dy:d);
%
%% Build base initial conditions once (reset each CFL run from these)
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
% Current in m/day (converted once, reset each run)
current_base = zeros(n+1, m+1, 2);
current_base(:,:,1) = mask .* u_max .* (1 - (2.*Y/(2*d)).^2);
current_base(:,:,2) = mask .* v_max;
current_base = current_base .* day;   % m/s -> m/day

%% =========================================================
%  CFL loop
% ==========================================================
for ci = 1:numel(cfl_targets)

    target_cfl = cfl_targets(ci);

    %% Time discretization derived from target CFL
    %  CFL = fish_speed_ms * dt_s / dx_m  =>  dt_s = CFL * dx_m / fish_speed_ms
    fish_speed = fish_speed_ms * day;           % m/day  (local, reset each run)
    dt_s       = target_cfl * dx_m / fish_speed_ms;   % seconds
    dt         = dt_s / day;                    % days
    num_hours  = dt * 24;
    iters      = ceil(days * 24 / num_hours);
    drawevery  = max(1, floor(num_hours));
    food_consumption_dt = eating_factor * dt;

    fprintf('\n======= CFL = %.2f  |  dt = %.4f days (%.4f h)  |  iters = %d =======\n', ...
            target_cfl, dt, num_hours, iters);
    fprintf('Fish speed: %.1f m/day  |  dx = %.1f m\n', fish_speed, dx_m);

    %% Reset fields for this run
    Fish    = Fish_base;
    Food    = Food_base;
    current = current_base;

    current(:,:,1) = squeeze(current(:,:,1)) .* mask;
    current(:,:,2) = squeeze(current(:,:,2)) .* mask;
    Food = Food .* mask;
    Fish = Fish .* mask;

    %% Initial True Solution
    delta_x  = x0_fish + (0 * dt_s * fish_speed_ms / 1000);
    True_Fish = makeFood('humpfish', Fish_base, X, Y, delta_x, 0, 1, 0);
    True_Fish(True_Fish < epsilon) = 0;

    apparent_current = current;

    %% Fish mass tracking: mass_total(k) [g] at time mass_time(k) [days]
    mass_time     = (0:iters) * dt;
    mass_total    = zeros(1, iters+1);
    mass_total(1) = sum(Fish, 'all') * vol_m;

    %% Iteration loop
    for it = 1:iters
        Fish = Fish .* mask;

        Fish_init = Fish;
        Food_init = Food;

        [Fish, Food, apparent_current_full] = ...
                      RandomSemiLagrangianAdvectFish( Fish, ...
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

        if do_smooth == 1
            Fish = smooth2D(Fish, 2);
        end

        fprintf('Iteration: %5d/%5d\n', it, iters);

        % Advance true solution
        delta_x  = x0_fish + (it * (fish_speed * dt) / 1000);
        True_Fish = makeFood('humpfish', Fish_init, X, Y, delta_x, 0, 1, 0);
        True_Fish(True_Fish < epsilon) = 0;

        % Record total fish mass after this step
        mass_total(it+1) = sum(Fish, 'all') * vol_m;
    end

    %% Error Calculations
    F_vec = Fish(:);
    A_vec = True_Fish(:);

    Error_Total = 1/(m*n) * sum((F_vec - A_vec).^2);
    corrFA = corr(F_vec, A_vec);
    stdF   = std(F_vec);    stdA  = std(A_vec);
    meanF  = mean(F_vec);   meanA = mean(A_vec);
    E_diss = (stdF - stdA)^2 + (meanF - meanA)^2;
    E_disp = 2*(1 - corrFA)*stdF*stdA;

    fprintf('CFL=%.2f:  MSE=%.6f  Corr=%.4f  E_diss=%.6f  E_disp=%.6f\n', ...
            target_cfl, Error_Total, corrFA, E_diss, E_disp);

    mass_rel_change = (mass_total(end) - mass_total(1)) / mass_total(1);
    fprintf('Fish mass: M0 = %.6e g, M_end = %.6e g, relative change = %+.3e\n', ...
            mass_total(1), mass_total(end), mass_rel_change);

    %% Save
    cfl_tag      = sprintf('cfl%03d', round(target_cfl * 100));
    savefilename = ['test1_Linear_RA_', cfl_tag, '.mat'];
    save(savefilename, 'Fish', 'X', 'Y', 'True_Fish', 'dx_m', 'dy_m', 'dt', 'dt_s', 'target_cfl', ...
         'mass_time', 'mass_total');
    fprintf('Saved: %s\n', savefilename);

end % CFL loop
