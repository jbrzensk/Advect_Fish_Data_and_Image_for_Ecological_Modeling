function[ Predator, Prey, apparent_current_full ] = SemiLagrangianAdvectFish( Predator,...
                                                                    Prey, ...
                                                                    current, ...
                                                                    food_consumption_dt, ...
                                                                    dt, ...
                                                                    dx_m, ...
                                                                    dy_m, ...
                                                                    fish_speed, ...
                                                                    mask, ...
                                                                    m, ...
                                                                    n)
%SEMILAGRANGIANADVECTFISH One time step of semi-Lagrangian fish advection.
%
%   [Predator, Prey, apparent_current_full] = ...
%       SemiLagrangianAdvectFish(Predator, Prey, current, food_consumption_dt, ...
%                                dt, dx_m, dy_m, fish_speed, mask, m, n)
%
%   Steps:
%     1. The fish eat food (eatFoodDomain), and the food is advected by the
%        current (AdvectFood) and masked.
%     2. The fractional increase in food toward each neighbour
%        (percentMoreFoodFourWay) and the swim speeds corrected for the
%        current (ApparentCurrentFull) are computed.
%     3. Every cell with fish sends a fraction |speed|*dt/dx of its fish to
%        each neighbour with more food (moveFishFunction). This fraction is
%        the local CFL number, and the total is capped at the fish in the
%        cell.
%     4. The fluxes from all cells are added to the fish field, so fish
%        biomass is conserved.
%
% INPUTS:
%   Predator            : [n+1 x m+1] fish concentration [g/m^3]
%   Prey                : [n+1 x m+1] food concentration [g/m^3]
%   current             : [n+1 x m+1 x 2] current, (:,:,1) = u in x,
%                         (:,:,2) = v in y [m/day]
%   food_consumption_dt : fraction of the food eaten per unit of fish in
%                         one step (passed to eatFoodDomain)
%   dt                  : time step [days]
%   dx_m, dy_m          : grid spacing in x and y [m]
%   fish_speed          : maximum swim speed of the fish [m/day]
%   mask                : [n+1 x m+1] 1 for ocean, 0 for land
%   m, n                : number of grid intervals in x and y
%
% OUTPUTS:
%   Predator              : [n+1 x m+1] fish concentration after the step
%   Prey                  : [n+1 x m+1] food concentration after the step
%   apparent_current_full : [n+1 x m+1 x 4] swim speed toward each
%                           neighbour (Up, Down, Right, Left) [m/day]
%
% NOTES:
%   - Any time unit works if current, fish_speed and dt use the same one.
%
% Brzenski, Jared
% 06/30/2025
% -------------------------------------------------------------------------
%% The core of the advection fish code.
    % COBALT Inputs
    % Eat the food
    Prey = eatFoodDomain( Prey, Predator, food_consumption_dt );
    % advect the food
    Prey = AdvectFood( Prey, current, dt, dx_m, dy_m );
    % mask off any advected food
    Prey = Prey .* mask;  
    %
    % Calculate the percent change between cells
    percent_more_food_full = percentMoreFoodFourWay( Prey );
    % Calculate the apparent current, taking swim speed and current
    apparent_current_full = ApparentCurrentFull( current, fish_speed, percent_more_food_full);
    %
    Fish_OG = Predator; % Hold original fish count
    %
    Flux = zeros(size(Predator));
    %
    for j=1:m+1
        for i=1:n+1
            idx = [i j];
            if Fish_OG(i,j) > 0
                
               Flux =  moveFishFunction( Predator, Flux, Prey, apparent_current_full, ...
                                          idx, ...
                                          dt, ...
                                          dx_m, ...
                                          dy_m );
            end
        end
    end
    Predator = Predator + Flux;
    %
end