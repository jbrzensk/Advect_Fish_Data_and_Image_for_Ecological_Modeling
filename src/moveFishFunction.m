function [ Flux ] = moveFishFunction( Predator, Flux, Prey, apparent_current_full, idx, dt, dx_m, dy_m )
%MOVEFISHFUNCTION Fish flux out of one cell toward neighbours with more food.
%
%   Flux = moveFishFunction(Predator, Flux, Prey, apparent_current_full, idx, dt, dx_m, dy_m)
%
%   Finds the neighbours of cell idx with more food (MoreFoodFourWay) and
%   passes them, with the swim speeds of the cell, to semiLagrangianFish,
%   which adds the fish moving out of cell idx to Flux.
%
% INPUTS:
%   Predator              : [ny x nx] fish concentration [g/m^3]
%   Flux                  : [ny x nx] change in fish concentration from
%                           movement this step, accumulated over cells
%   Prey                  : [ny x nx] food concentration [g/m^3]
%   apparent_current_full : [ny x nx x 4] swim speed toward each neighbour
%                           (Up, Down, Right, Left) [m/day]
%   idx                   : [i, j] row and column of the cell
%   dt                    : time step [days]
%   dx_m, dy_m            : grid spacing in x and y [m]
%
% OUTPUT:
%   Flux                  : Flux with the movement out of cell idx added
%
% Brzenski, Jared
% 06/30/2025
% -------------------------------------------------------------------------
%% The core of the move fish code.
i = idx(1);
j = idx(2);

[food_dir] = MoreFoodFourWay( Prey, idx );

apparent_current_single = reshape(apparent_current_full(i,j,:), [4, 1]);

Flux = semiLagrangianFish( Predator, Flux, [i,j], food_dir, ...
                         apparent_current_single, ...
                         dt, dx_m, dy_m);
end