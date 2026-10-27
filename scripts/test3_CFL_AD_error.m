%% CFL Test, linear advection+diffusion of fish
%  Fish follow food
%
% Test Case 1: Linear advection for 3 days.
% No consumption, no speed adjustment, just
% if Conc_n > Conc_m, swim towards Cell_n at max_speed
%
% The original concentration should be advected byL
% 3days * (1.0m/s * 60*60*24) meters
%
% DIFFUSION: Average diffusion value of 1000 m^2/s
%
% BY: Jared Brzenski
%
% Units:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This first section has the editable quantities:
% iters            : number of steps to run
% drawevery        : how often output should be plotted.
% eating_factor    : what percent of the fish body mass (concentration)
%                    is eaten from the food source in a time step.
% do_smooth        : artificial diffusion. A moving average of X cells
%                    which conserves the quantities smoothed.
%
clc
close all
clearvars

addpath('../src')

savefilename = 'CFL_AD_error.mat';

kappa = 1000; % in m^2/s

speed_tries = linspace(0.05,1,20);
time_tries = linspace(1,12,47);

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
parfor ii=1:speeds
    for jj=1:times
        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n' );
        fprintf('Interation: %4.0f out of %5.0f \n',jj+(ii-1)*speeds, speeds*times );
        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n' );
        fish_speed = speed_tries(ii); 
        num_hours = time_tries(jj);
        
        %% Time / draw
        % days = 200;
        days = 3;      % sim time
        %num_hours = 2; % dt
        
        iters = ceil(days * 24 / num_hours);             
        drawevery = floor(num_hours);           % update plot after this many iterations
        %
        %% System Parameters
        % Fish Swimming Speed
        % Here BATS Tuna Speed 86-260km/day ~ 1m/s to 3m/s
        % Biomass is approximately 1gww/m^3 ~ 
        % Tuna consume ~5% of their mass per day
        %fish_speed = 1.0;       % speed in meters/second
        eating_factor = 0.00;   % percentage of food eaten by fish vs concentration
        %
        fish_init = 1.0;        % fish biomass in g / m^3
        food_init = 1000;       % food biomass in g / m^3
        %
        % Currents
        u_max = 0.;           % current speed in meters/second
        v_max = 0.0;
        % Time Discretization
        dt   = num_hours/24;    % dt in days
        dt_s = dt*24*60*60;     % dt in seconds
        day  = 24*60*60;        % seconds in a day ( for dimensionalizing )
        %
        % diffusion kappa dimension to 
        kappa_dt = kappa * dt_s;
        %
        epsilon = 10e-15;      % tiny value, 
        %
        do_smooth = 0;           % do any smoothing
        %
        %% Food scenario, gradient, blob
        x0 = 800;
        y0 = 0;
        amp = food_init;
        other_conc = 0.1*food_init;
        %
        %% Grid setup
        % Number of CELLS in each direction
        m = 100;
        n = 50;
        % Domain's start and end locations ( distances )
        a = 0;
        b = 1000;
        c = -250.0;
        d = 250.0;
        % Spatial step sizes
        dx = (b-a)/m;         % kilometers
        dy = (d-c)/n;         % kilometers
        dx_m = dx * 1000.0;   % dx in meters
        dy_m = dy * 1000.0;   % dy in meters
        %
        % Volumes
        vol = dx * dy * 0.1;             % volume of water (km^3), assuming 100m thick
        vol_m = vol * ( 1000^2 * 100 );  % volume  of water (meters^3)
        %
        length_u = (m+1)*n;
        length_v = m*(n+1);
        length_c = (m+2)*(n+2);
        %
        % Grid
        [X, Y] = meshgrid( a : dx : b, c : dy : d );
        %
        %% Initial Density Generation for Fish and Food
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameters
        Fish = fish_init * ones(n+1, m+1);
        Fish = zeroCornerMatrix(Fish);
        Fish_init = Fish;

        % Fish using Fodd hump!
        x0_fish = 300; % km
        Fish = makeFood('humpfish', Fish_init, X, Y, x0_fish, 0, 1, 0);
        
        Fish(Fish < epsilon) = 0;
        
        %% True Solution after 3 days in km
        delta_x = x0_fish + (0*dt_s*fish_speed/1000);
        True_Fish = makeFood('humpfish', Fish_init, X, Y, delta_x, 0, 1, 0);
        True_Fish(True_Fish < epsilon) = 0;
        
        %% Initial Food Blob Generation
        %
        Food = food_init * ones(n+1, m+1);
        Food = makeFood('gradient', Food, X, Y, x0, y0, amp, other_conc);
        Food = fliplr(Food); % Flip gadient left to right
        Food = Food .* food_init;
        %
        %% Island Maskingcobalt%hp_ingest_nlgz
        mask = ones(size(Food));
        
        Fish = mask .* Fish;
        Food = mask .* Food; 
        %% Initial Current Field
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        current = zeros(n+1,m+1,2);
        % U init
        current(:,:,1) = mask .* u_max;
        % Parabolci inlet velocity
        current(:,:,1) = mask .* u_max .*( 1 - (2.*Y/(2*d)).^2 );
        % V init
        current(:,:,2) = mask .* v_max;
        %
        %% Convert to day from seconds
        % Convert to m per day
        current = current .* day;
        % Fish speed to m per day
        fish_speed = fish_speed .* day;
        % Consumption per time step
        food_consumption_dt = eating_factor * dt;
        fprintf('Max Fish Swim Speed: %4.5f m/dt\n', fish_speed*dt);
        
        %
        fprintf('Max U Current: %4.5f m/day\n', max(current(:)) );
        fprintf('Max U Current: %4.5f m/dt\n', max(current(:))*dt );
        fprintf('Max Fish Swim Speed: %4.5f m/day\n', fish_speed);
        fprintf('Max Fish Swim Speed: %4.5f m/dt\n', fish_speed*dt);
        fprintf('dx: %4.5f [m]\n', dx_m);
        fprintf('dy: %4.5f [m]\n', dy_m);
        %
        
        %% Mask Values
        current(:,:,1) = squeeze(current(:,:,1)) .* mask;
        current(:,:,2) = squeeze(current(:,:,2)) .* mask;
        
        apparent_current = current;
        
        Food = Food .* mask;
        Fish = Fish .* mask;
        
        %
        %% Fish mass tracking
        t_hist    = (0:iters) * dt;
        m_hist    = zeros(1, iters+1);
        m_hist(1) = sum(Fish, 'all') * vol_m;

        %% Iteration Starts
        for it = 1 : iters
            % Mask fo fish
            Fish = Fish .* mask;
        
            % Conserve Fish and Food quantity
            Fish_init = Fish;
            Food_init = Food;
            %
            [ Fish, Food, apparent_current_full ] = SemiLagrangianAdvectFish( Fish,...
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
        
            Fish = laplacian2D( Fish, dx_m, dy_m, kappa_dt );


            delta_x = x0_fish + (it*(fish_speed*dt)/1000);
            True_Fish = makeFood('humpfish', Fish_init, X, Y, delta_x, 0, 1, 0);
            True_Fish(True_Fish < epsilon) = 0;

            % Record total fish mass after this step
            m_hist(it+1) = sum(Fish, 'all') * vol_m;
            %
        end
        
        mass_time{ii,jj}  = t_hist;
        mass_total{ii,jj} = m_hist;
        fprintf('Fish mass relative change: %+.3e\n', (m_hist(end) - m_hist(1)) / m_hist(1));

        %% Error Calculations
        cfl(ii,jj) = fish_speed * dt / dx_m;
        fprintf('CFL: %3.5e \n',cfl(ii,jj) );

        F_vec = Fish(:);
        A_vec = True_Fish(:);
        
        E_total(ii,jj) = 1/(m*n) * sum( (F_vec(:) - A_vec(:)).^2 );
        
        corrFA = corr(F_vec, A_vec);
        
        stdF = std(F_vec);
        stdA = std(A_vec);
        
        meanF = mean(F_vec);
        meanA = mean(A_vec);
        
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

sorted_vals     = sort(scaled_E_total(:));
thresh_idx      = round(0.25 * numel(sorted_vals));
threshold_value = sorted_vals(thresh_idx);

hold on;
contour(X_plt, Y_plt, scaled_E_total', [threshold_value threshold_value], ...
        'r', 'LineWidth', 3, 'DisplayName', 'min 25\% region');
contour_levels = linspace(0.1, 2, 20);
[C5, h5] = contour( X_plt, Y_plt, cfl', contour_levels, 'k', 'LineWidth', 1,...
                   'DisplayName','cfl');
clabel(C5, h5, 'FontSize', clabelfont, 'Color', 'k', 'LabelSpacing', 400);
hold off;

pbaspect([1 1 1]);
xlabel('$S_{\max}$ [m/s]', 'Interpreter', 'latex', 'FontSize', labelfont);
ylabel('$\Delta t$ [hr]',  'Interpreter', 'latex', 'FontSize', labelfont);
title('Advection + Diffusion: Total Error with CFL Contours', ...
      'Interpreter', 'latex', 'FontSize', titlefont);
set(gca, 'FontSize', labelfont, 'TickLabelInterpreter', 'latex');

save(savefilename, 'X_plt','Y_plt','scaled_E_total','cfl', 'scaled_E_diss', 'scaled_E_disp', ...
     'mass_time', 'mass_total', 'mass_rel_change')