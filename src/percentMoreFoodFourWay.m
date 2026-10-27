function [ percent_more_food_full ] = percentMoreFoodFourWay( Food )
%PERCENTMOREFOODFOURWAY Fractional increase in food toward each neighbour.
%
%   percent_more_food_full = percentMoreFoodFourWay(Food)
%
%   For every cell and each of its four neighbours, returns
%   min(F_neighbour / F_cell, 2) - 1 when the neighbour has more food, and
%   0 otherwise. The value is between 0 and 1, i.e. capped at a 100%
%   increase. Neighbours outside the grid give 0.
%
% INPUTS:
%   Food                   : [ny x nx] food concentration
%
% OUTPUT:
%   percent_more_food_full : [ny x nx x 4] fractional increase toward
%                            Up (row-1), Down (row+1), Right (col+1),
%                            Left (col-1)
%
% Brzenski, Jared
% 06/30/2025
% -------------------------------------------------------------------------
%% Core percent more food code
    % Define neighbors relative to (i, j)
    neighbors = [-1, 0;  % Up
                  1, 0;  % Down
                  0, 1;  % Right
                  0, -1];% Left

    [n, m] = size(Food);

    percent_more_food_full = zeros(n,m,4);

    percent_more_food = [ 0 0 0 0 ];

    for j=1:m
        for i=1:n
            current_val = Food(i,j);
            % Reset food vector each loop
            percent_more_food = [ 0 0 0 0 ];
            % Loop through each direction ( find more food )
            for k = 1:4
                ni = i + neighbors(k, 1); % Neighbor row
                nj = j + neighbors(k, 2); % Neighbor column
                
                % Check bounds and compare food values
                if ni >= 1 && ni <= n && nj >= 1 && nj <= m % Valid neighbor
                    if Food(ni, nj) > current_val
                        percent_more_food(k) = min( Food(ni,nj)/current_val, 2 ) - 1;
                    end
                end
            end
            percent_more_food_full(i,j,:) = percent_more_food;
        end
    end
end
