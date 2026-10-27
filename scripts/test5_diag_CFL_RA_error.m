%% CFL Test, diagonal random advection of fish (no diffusion)
%  Fish follow a diagonal food gradient (x+y direction) via random
%  semi-Lagrangian advection
%
% Test Case 5RA: Diagonal random advection CFL error sweep.
% No consumption, no speed adjustment, no diffusion.
%
% Diagonal IC: fish start at (500, 0) on a 200x200 grid spanning
% [0,1000] x [-500,500] km.  Food gradient is diagonal (Y+250+X)/1000,
% so fish swim at 45 degrees. True solution advances at fish_speed/sqrt(2)
% in both x and y.
%
% Tests the same CFL sweep as test3_CFL_RA_error.m but for diagonal advection.
%
% Sweeps over fish_speed and dt to compute dissipation, dispersion,
% and total error relative to the true solution at 2.5 days.
%
% BY: Jared Brzenski
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clearvars

addpath('../src')

% Random number set
rng(42, 'twister');   % 'twister' = Mersenne Twister, MATLAB's default algorithm

savefilename = 'CFL_Diag_RA_error.mat';

speed_tries = linspace(0.05,1,20);
time_tries  = linspace(1,12,47);

[ X, Y ] = meshgrid( speed_tries, time_tries );

times   = length(time_tries);
speeds  = length(speed_tries);
E_total = zeros(speeds, times);
E_diss  = zeros(speeds, times);
E_disp  = zeros(speeds, times);
cfl     = zeros(speeds, times);
% Fish mass history per run (iters differs per run, so cell arrays):
% mass_total{ii,jj}(k) [g] at time mass_time{ii,jj}(k) [days]
mass_time  = cell(speeds, times);
mass_total = cell(speeds, times);

tic
for ii=1:speeds
 parfor jj=1:times
        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n' );
        fprintf('Interation: %4.0f out of %5.0f \n',jj+(ii-1)*speeds, speeds*times );
        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n' );
        fish_speed = speed_tries(ii);
        num_hours  = time_tries(jj);

        %% Time / draw
        days = 2.5;       % sim time
        iters     = ceil(days * 24 / num_hours);
        drawevery = floor(num_hours);

        %% System Parameters
        eating_factor = 0.00;
        fish_init     = 1.0;
        food_init     = 1000;
        u_max = 0.;
        v_max = 0.0;
        dt    = num_hours/24;   % dt in days
        dt_s  = dt*24*60*60;    % dt in seconds
        day   = 24*60*60;
        epsilon = 10e-15;
        do_smooth = 0;

        %% Grid setup
        m = 100;
        n = 100;
        a = 0;     b = 1000;
        c = -500.0; d = 500.0;
        dx   = (b-a)/m;
        dy   = (d-c)/n;
        dx_m = dx * 1000.0;
        dy_m = dy * 1000.0;

        vol   = dx * dy * 0.1;
        vol_m = vol * ( 1000^2 * 100 );

        length_u = (m+1)*n;
        length_v = m*(n+1);
        length_c = (m+2)*(n+2);

        [X_grid, Y_grid] = meshgrid( a : dx : b, c : dy : d );

        %% Initial Conditions
        Fish = fish_init * ones(n+1, m+1);
        Fish = zeroCornerMatrix(Fish);
        Fish_init = Fish;

        x0_fish = 500;
        y0_fish = 0;
        Fish = makeFood('humpfish', Fish_init, X_grid, Y_grid, x0_fish, y0_fish, 1, 0);
        Fish(Fish < epsilon) = 0;

        %% True Solution (no advection at t=0)
        True_Fish = makeFood('humpfish', Fish_init, X_grid, Y_grid, x0_fish, y0_fish, 1, 0);
        True_Fish(True_Fish < epsilon) = 0;

        %% Initial Food (diagonal gradient, no makeFood)
        Food = ((Y_grid + 250) + X_grid) / 1000;
        Food = Food .* food_init;

        %% Mask
        mask = ones(size(Food));
        Fish = mask .* Fish;
        Food = mask .* Food;

        %% Current Field
        current = zeros(n+1,m+1,2);
        current(:,:,1) = mask .* u_max;
        current(:,:,1) = mask .* u_max .* ( 1 - (2.*Y_grid/(2*d)).^2 );
        current(:,:,2) = mask .* v_max;

        %% Unit conversions
        current    = current .* day;
        fish_speed = fish_speed .* day;
        food_consumption_dt = eating_factor * dt;

        current(:,:,1) = squeeze(current(:,:,1)) .* mask;
        current(:,:,2) = squeeze(current(:,:,2)) .* mask;

        apparent_current = current;
        Food = Food .* mask;
        Fish = Fish .* mask;

        %% Fish mass tracking
        t_hist    = (0:iters) * dt;
        m_hist    = zeros(1, iters+1);
        m_hist(1) = sum(Fish, 'all') * vol_m;

        %% Iteration loop
        root2 = sqrt(2);
        for it = 1 : iters
            Fish = Fish .* mask;

            Fish_init = Fish;
            Food_init = Food;

            [ Fish, Food, apparent_current_full ] = RandomSemiLagrangianAdvectFish( Fish, ...
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

            % No diffusion in this test

            distance  = it*(fish_speed*dt)/1000/root2;
            delta_x   = x0_fish + distance;
            delta_y   = y0_fish + distance;
            True_Fish = makeFood('humpfish', Fish_init, X_grid, Y_grid, delta_x, delta_y, 1, 0);
            True_Fish(True_Fish < epsilon) = 0;

            % Record total fish mass after this step
            m_hist(it+1) = sum(Fish, 'all') * vol_m;
        end

        mass_time{ii,jj}  = t_hist;
        mass_total{ii,jj} = m_hist;
        fprintf('Fish mass relative change: %+.3e\n', (m_hist(end) - m_hist(1)) / m_hist(1));

        %% Error Calculations
        cfl(ii,jj) = fish_speed * dt / dx_m;
        fprintf('CFL: %3.5e \n', cfl(ii,jj));

        F_vec = Fish(:);
        A_vec = True_Fish(:);

        E_total(ii,jj) = 1/(m*n) * sum( (F_vec(:) - A_vec(:)).^2 );

        corrFA = corr(F_vec, A_vec);
        stdF   = std(F_vec);   stdA  = std(A_vec);
        meanF  = mean(F_vec);  meanA = mean(A_vec);

        E_diss(ii,jj) = ( stdF - stdA )^2 + ( meanF - meanA )^2;
        E_disp(ii,jj) = 2*(1-corrFA)*stdF*stdA;
    end
end
toc
%% Mass Conservation Summary
mass_rel_change = cellfun(@(M) (M(end) - M(1)) / M(1), mass_total);
fprintf('Max |relative fish mass change| over all runs: %.3e\n', ...
        max(abs(mass_rel_change(:))));

%% Scale for plotting
scaled_E_diss  = E_diss  .* 10000;
scaled_E_disp  = E_disp  .* 1000;
scaled_E_total = E_total .* 1000;

%% Plotting — publication quality
labelfont    = 14;
titlefont    = 15;
clabelfont   = 11;
lineweight   = 1.5;
num_contours = 10;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');

[ X_plt, Y_plt ] = meshgrid( speed_tries, time_tries );

%% Figure 1: 2×2 panel — error components and CFL
fig1 = figure('Units', 'inches', 'Position', [1 1 8 7]);
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'tight');

% (a) Dissipation
ax1 = nexttile;
[C1, h1] = contour( X_plt, Y_plt, scaled_E_diss', num_contours, 'LineWidth', lineweight );
clabel(C1, h1, 'FontSize', clabelfont, 'Color', 'k', 'LabelSpacing', 400);
axis tight; pbaspect([1 1 1]);
xlabel('$S_{\max}$ [m/s]', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\Delta t$ [hr]',  'Interpreter', 'latex', 'FontSize', labelfont);
title('(a) Dissipation $\epsilon_{diss}$ ($\times 10^{-4}$)', ...
      'Interpreter', 'latex', 'FontSize', titlefont);
set(ax1, 'FontSize', labelfont, 'TickLabelInterpreter', 'latex');

% (b) Dispersion
ax2 = nexttile;
[C2, h2] = contour( X_plt, Y_plt, scaled_E_disp', num_contours, 'LineWidth', lineweight );
clabel(C2, h2, 'FontSize', clabelfont, 'Color', 'k', 'LabelSpacing', 400);
axis tight; pbaspect([1 1 1]);
xlabel('$S_{\max}$ [m/s]', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\Delta t$ [hr]',  'Interpreter', 'latex', 'FontSize', labelfont);
title('(b) Dispersion $\epsilon_{disp}$ ($\times 10^{-3}$)', ...
      'Interpreter', 'latex', 'FontSize', titlefont);
set(ax2, 'FontSize', labelfont, 'TickLabelInterpreter', 'latex');

% (c) Total Error
ax3 = nexttile;
[C3, h3] = contour( X_plt, Y_plt, scaled_E_total', num_contours, 'LineWidth', lineweight );
clabel(C3, h3, 'FontSize', clabelfont, 'Color', 'k', 'LabelSpacing', 400);
axis tight; pbaspect([1 1 1]);
xlabel('$S_{\max}$ [m/s]', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\Delta t$ [hr]',  'Interpreter', 'latex', 'FontSize', labelfont);
title('(c) Total Error $\epsilon_{tot}$ ($\times 10^{-3}$)', ...
      'Interpreter', 'latex', 'FontSize', titlefont);
set(ax3, 'FontSize', labelfont, 'TickLabelInterpreter', 'latex');

% (d) CFL Number
ax4 = nexttile;
contour_levels = linspace(0.1, 2, 20);
[C4, h4] = contour( X_plt, Y_plt, cfl', contour_levels, 'LineWidth', lineweight );
clabel(C4, h4, 'FontSize', clabelfont, 'Color', 'k', 'LabelSpacing', 400);
axis tight; pbaspect([1 1 1]);
xlabel('$S_{\max}$ [m/s]', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\Delta t$ [hr]',  'Interpreter', 'latex', 'FontSize', labelfont);
title('(d) CFL Number', 'Interpreter', 'latex', 'FontSize', titlefont);
set(ax4, 'FontSize', labelfont, 'TickLabelInterpreter', 'latex');

%% Figure 2: Total error (filled) with CFL overlay
fig2 = figure('Units', 'inches', 'Position', [1 1 8 8]);
contourf(X_plt, Y_plt, scaled_E_total', 10, 'LineStyle', 'none', 'DisplayName', 'error');
colormap(fig2, parula);
cb = colorbar;
ylabel(cb, '$\epsilon_{tot}\ (\times 10^{-3})$', ...
       'Interpreter', 'latex', 'FontSize', labelfont);
set(cb, 'TickLabelInterpreter', 'latex', 'FontSize', labelfont-1);

% 25 percent region
sorted_vals     = sort(scaled_E_total(:));
thresh_idx      = round(0.25 * numel(sorted_vals));
threshold_value = sorted_vals(thresh_idx);
hold on;
contour(X_plt, Y_plt, scaled_E_total', [threshold_value threshold_value], ...
        'r', 'LineWidth', 3, 'DisplayName', 'min 25\% region');

% CFL values
contour_levels = linspace(0.1, 2, 20);
[C5, h5] = contour( X_plt, Y_plt, cfl', contour_levels, 'k', ...
                   'LineWidth', 1, 'DisplayName', 'cfl' );
clabel(C5, h5, 'FontSize', clabelfont, 'Color', 'k', 'LabelSpacing',400);
hold off;

pbaspect([1 1 1]);
xlabel('$S_{\max}$ [m/s]', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\Delta t$ [hr]',  'Interpreter', 'latex', 'FontSize', labelfont);
title('Diagonal Random Advection: Total Error with CFL Contours', ...
      'Interpreter', 'latex', 'FontSize', titlefont);
set(gca, 'FontSize', labelfont, 'TickLabelInterpreter', 'latex');
legend

save(savefilename, 'X_plt','Y_plt','scaled_E_total','cfl', 'scaled_E_diss', 'scaled_E_disp', ...
     'mass_time', 'mass_total', 'mass_rel_change')
