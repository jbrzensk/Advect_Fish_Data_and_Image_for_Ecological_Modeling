function [Food_base, Fish_base, True_Fish_final, X, Y] = horizontal_initial_conditions()

%% Fixed parameters (same for every CFL run)
days          = 3;      % sim time [days]
fish_speed_ms = 1.0;    % fish speed [m/s]  — converted to m/day inside loop
fish_init     = 1.0;    % fish biomass [g/m^3]
food_init     = 1000;   % food biomass [g/m^3]
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
%% Grid (fixed for all runs)
m  = 100;   n  = 50;
a  = 0;     b  = 1000;
c  = -250.0; d = 250.0;
dx   = (b-a)/m;           % km
dy   = (d-c)/n;           % km
%
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
delta_x_final   = x0_fish + (fish_speed_ms * day * days) / 1000;   % km
True_Fish_final = makeFood('humpfish', Fish_base, X, Y, delta_x_final, 0, 1, 0);
True_Fish_final(True_Fish_final < epsilon) = 0;