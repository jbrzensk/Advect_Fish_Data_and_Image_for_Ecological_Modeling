function[ Fish, Food, apparent_current_full, percent_more_food_full ] = RandomSemiLagrangianAdvectFish( Fish,...
                                                                                                          Food, ...
                                                                                                          current, ...
                                                                                                          food_consumption_dt, ...
                                                                                                          dt, ...
                                                                                                          dx_m, ...
                                                                                                          dy_m, ...
                                                                                                          fish_speed, ...
                                                                                                          mask, ...
                                                                                                          m, ...
                                                                                                          n)
%RANDOMSEMILAGRANGIANADVECTFISH One time step of random semi-Lagrangian fish advection (RA).
%
%   [Fish, Food, apparent_current_full, percent_more_food_full] = ...
%       RandomSemiLagrangianAdvectFish(Fish, Food, current, food_consumption_dt, ...
%                                      dt, dx_m, dy_m, fish_speed, mask, m, n)
%
%   Same steps as SemiLagrangianAdvectFish, with random swim directions:
%     1. The fish eat food (eatFoodDomain), and the food is advected by the
%        current (AdvectFood) and masked.
%     2. The fractional increase in food toward each neighbour
%        (percentMoreFoodFourWay) and the swim speeds corrected for the
%        current (ApparentCurrentFull) are computed.
%     3. Random direction weights are drawn for every cell (randomSpeed).
%     4. Every cell with fish sends fish to its neighbours with more food
%        (RandomMoveFishFunction), with the speeds scaled by the weights.
%
% INPUTS:
%   Fish                : [n+1 x m+1] fish concentration [g/m^3]
%   Food                : [n+1 x m+1] food concentration [g/m^3]
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
%   Fish                   : [n+1 x m+1] fish concentration after the step
%   Food                   : [n+1 x m+1] food concentration after the step
%   apparent_current_full  : [n+1 x m+1 x 4] swim speed toward each
%                            neighbour (Up, Down, Right, Left) [m/day]
%   percent_more_food_full : [n+1 x m+1 x 4] fractional increase in food
%                            toward each neighbour
%
% NOTES:
%   - Any time unit works if current, fish_speed and dt use the same one.
%   - Random numbers come from MATLAB's global stream (through
%     circ_vmrnd), so call rng first to make a run repeatable.
%
% Brzenski, Jared
% 03/09/2026
% -------------------------------------------------------------------------
%% The core of the advection fish code.
    % COBALT Inputs
    % Eat the food
    Food = eatFoodDomain( Food, Fish, food_consumption_dt );
    % advect the food
    Food = AdvectFood( Food, current, dt, dx_m, dy_m );
    % mask off any advected food
    Food = Food .* mask;  
    %
    % Calculate the percent change between cells
    percent_more_food_full = percentMoreFoodFourWay( Food );
    %percent_more_food_full = percentMoreFoodFourWaySigned( Food );
    % Calculate the apparent current, taking swim speed and current
    apparent_current_full = ApparentCurrentFull( current, fish_speed, percent_more_food_full);
    %
    rand_weighting = randomSpeed( apparent_current_full, fish_speed, percent_more_food_full );
    %
    %save('Test_arrays_for_random.mat', 'apparent_current_full','Fish','Food','percent_more_food_full');
    %
    Fish_OG = Fish; % Hold original fish count
    %
    Flux = zeros(size(Fish));
    %
    for j=1:m+1
        for i=1:n+1
            idx = [i j];
            if Fish_OG(i,j) > 0
                
               Flux =  RandomMoveFishFunction( Fish, Flux, Food, apparent_current_full, ...
                                              rand_weighting, ...
                                              idx, ...
                                              dt, ...
                                              dx_m, ...
                                              dy_m );
            end
        end
    end
    Fish = Fish + Flux;
    %
end