function [ dir ] = MoreFoodFourWay( Food, idx )
%MOREFOODFOURWAY Find which neighbours of one cell have more food.
%
%   dir = MoreFoodFourWay(Food, idx)
%
%   Compares the food in cell idx with its four neighbours. Neighbours
%   outside the grid count as having no more food.
%
% INPUTS:
%   Food  : [ny x nx] food concentration
%   idx   : [i, j] row and column of the cell
%
% OUTPUT:
%   dir   : [1 x 4] 1 where the neighbour has more food than cell idx,
%           0 otherwise. Order: Up (row-1), Down (row+1), Right (col+1),
%           Left (col-1)
%
% Brzenski, Jared
% 06/30/2025
% -------------------------------------------------------------------------
%% Core percent more food code
% Finds the direction of more food ONLY! 
    dir = [ 0 0 0 0 ];
    
    percent_more_food = [ 0 0 0 0 ];

    i = idx(1);
    j = idx(2);

    current_val = Food(i,j);

    [n, m] = size(Food);
    
    % Define neighbors relative to (i, j)
    neighbors = [-1, 0;  % Up
                  1, 0;  % Down
                  0, 1;  % Right
                  0, -1];% Left

    % Loop through each direction ( find more food )
    for k = 1:4
        ni = i + neighbors(k, 1); % Neighbor row
        nj = j + neighbors(k, 2); % Neighbor column
        
        % Check bounds and compare food values
        if ni >= 1 && ni <= n && nj >= 1 && nj <= m % Valid neighbor
            if Food(ni, nj) > current_val
                dir(k) = 1;
            end
        end
    end
    
end
