function [ Flux ] = RandomMoveFishFunction( Fish, Flux, Food, apparent_current_full, rand_weighting, idx, dt, dx_m, dy_m )
%RANDOMMOVEFISHFUNCTION Random-advection fish flux out of one cell.
%
%   Flux = RandomMoveFishFunction(Fish, Flux, Food, apparent_current_full, rand_weighting, idx, dt, dx_m, dy_m)
%
%   Same as moveFishFunction, but also passes the random direction weights
%   of cell idx to RandomSemiLagrangianFish, which scales the swim speeds
%   by them.
%
% INPUTS:
%   Fish                  : [ny x nx] fish concentration [g/m^3]
%   Flux                  : [ny x nx] change in fish concentration from
%                           movement this step, accumulated over cells
%   Food                  : [ny x nx] food concentration [g/m^3]
%   apparent_current_full : [ny x nx x 4] swim speed toward each neighbour
%                           (Up, Down, Right, Left) [m/day]
%   rand_weighting        : [ny x nx x 4] random weights from randomSpeed
%   idx                   : [i, j] row and column of the cell
%   dt                    : time step [days]
%   dx_m, dy_m            : grid spacing in x and y [m]
%
% OUTPUT:
%   Flux                  : Flux with the movement out of cell idx added
%
% Brzenski, Jared
% 05/27/2025
% -------------------------------------------------------------------------
i = idx(1);
j = idx(2);

[food_dir] = MoreFoodFourWay( Food, idx );

rand_weight = squeeze( rand_weighting(i,j,:) );

apparent_current_single = squeeze(apparent_current_full(i,j,:));

Flux = RandomSemiLagrangianFish( Fish, Flux, [i,j], rand_weight, ...
                                 food_dir, apparent_current_single, ...
                                 dt, dx_m, dy_m);
end