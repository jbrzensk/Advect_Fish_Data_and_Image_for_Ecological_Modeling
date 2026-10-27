function [Flux] = semiLagrangianFish(conc_matrix, Flux, idx, dir, current, dt, dx_m, dy_m)
%SEMILAGRANGIANFISH Move fish out of one cell toward neighbours with more food.
%
%   Flux = semiLagrangianFish(conc_matrix, Flux, idx, dir, current, dt, dx_m, dy_m)
%
%   For each direction flagged in dir, a fraction |speed|*dt/distance of
%   the fish in cell idx moves to that neighbour, where distance is dy_m
%   for Up/Down and dx_m for Right/Left. This fraction is the local CFL
%   number. If the total would exceed the fish in the cell, all fractions
%   are scaled down so the cell empties exactly. The fish leaving are
%   subtracted from Flux at idx and added at the neighbours, so mass is
%   conserved.
%
% INPUTS:
%   conc_matrix : [ny x nx] fish concentration [g/m^3]
%   Flux        : [ny x nx] change in fish concentration from movement this
%                 step, accumulated over cells
%   idx         : [i, j] row and column of the cell
%   dir         : [1 x 4] 1 for neighbours with more food, 0 otherwise
%                 (Up, Down, Right, Left)
%   current     : [4 x 1] swim speed toward each neighbour [m/day]
%   dt          : time step [days]
%   dx_m, dy_m  : grid spacing in x and y [m]
%
% OUTPUT:
%   Flux        : Flux with the movement out of cell idx added
%
% NOTES:
%   - Fish move at most one cell per step, so for CFL > 1 they travel
%     slower than the swim speed.
%
% Brzenski, Jared
% 06/30/2025
% -------------------------------------------------------------------------
    % Directions matrix
    directions = [
                -1, 0;  % Up
                 1, 0;  % Down
                 0, 1;  % Right
                 0, -1];% Left  
    
    % Only do four core directions
    dir = dir(1:4);

    distance = [
                dy_m;
                dy_m;
                dx_m;
                dx_m;];

    [ n, m] = size( conc_matrix );

    % Current cell coordinates
    i = idx(1);
    j = idx(2);

    % Current cell concentration ( saves lookups later )
    cell_concentration = conc_matrix(i,j);

    % Calculate the number of active directions
    activeDirections = find( dir == 1);     % indices
    totalDirections = numel(activeDirections);  % number of indices

    % If no movement, return immediately
    if totalDirections == 0
        return;
    end
    
    % Speed swimming in each direction;
    speeds = dir' .* current;
    %
    proportion_out = abs(speeds) .* dt ./ distance;
    %
    conc_moving_out = proportion_out .* cell_concentration;
    %
    total_leaving = sum(conc_moving_out);
    %
    % If too much wants to leave, cap at cell capacity
    if total_leaving > cell_concentration
        % Scale amounts moving out to available concentration
        conc_moving_out = conc_moving_out * (cell_concentration / total_leaving );
        total_leaving = sum(conc_moving_out);
    end
    % OLD
    % conc_moving_out = min( proportion_out .* cell_concentration, cell_concentration/totalDirections );
    % Remove movers from the current cell
    %
    Flux(i,j) = Flux(i,j) - total_leaving; % update for fish leaving this timestep
    %
    %conc_matrix(i, j) = cell_concentration - total_leaving;

    % Update fish positions based on active directions
    for k = 1:totalDirections
        directionIndex = activeDirections(k); % Get the direction index
        ni = i + directions(directionIndex, 1); % Neighbor row
        nj = j + directions(directionIndex, 2); % Neighbor column

        % Ensure the neighbor is within bounds
        if ni >= 1 && ni <= n && nj >= 1 && nj <= m
            Flux(ni, nj) = Flux(ni, nj) + conc_moving_out(directionIndex);
            %conc_matrix(ni, nj) = conc_matrix_OG(ni, nj) + conc_moving_out(directionIndex);
        end
    end
end