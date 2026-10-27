function [ dir_swim_speed_full ] = ApparentCurrentFull( current, fish_speed, percent_more_food )
%APPARENTCURRENTFULL Swim speed in the four directions, corrected for the current.
%
%   dir_swim_speed_full = ApparentCurrentFull(current, fish_speed, percent_more_food)
%
%   In each cell, fish swim at fish_speed toward each neighbour that has
%   more food (percent_more_food > 0), and do not swim toward the others.
%   The current along each direction is added, and negative speeds are set
%   to zero. Where the current along a direction is at least fish_speed,
%   the speed in that direction is the excess, current - fish_speed.
%
% INPUTS:
%   current             : [ny x nx x 2] current, (:,:,1) = u in x,
%                         (:,:,2) = v in y
%   fish_speed          : scalar, maximum swim speed of the fish
%   percent_more_food   : [ny x nx x 4] fractional increase in food toward
%                         each neighbour (Up, Down, Right, Left), from
%                         percentMoreFoodFourWay
%
% OUTPUT:
%   dir_swim_speed_full : [ny x nx x 4] speed toward each neighbour
%                         (Up, Down, Right, Left), in the units of
%                         fish_speed
%
% NOTES:
%   - current and fish_speed must use the same units (m/day in the test
%     scripts).
%   - percent_more_food is only used to pick the directions with more
%     food. A speed that depends on the food gradient could be added here.
%
% Brzenski, Jared
% 06/30/2025
% -------------------------------------------------------------------------
%% Core Apparent Current Full Code
    % Directions matrix
    happy = 0;
    
    directions = [
                -1, 0;  % Up
                 1, 0;  % Down
                 0, 1;  % Right
                 0, -1];% Left

    [n, m, o] = size(current);

    dir_swim_speed_full = zeros(n,m,4);

    for i=1:n
        for j=1:m
   
            % Old
            % U = current(i,j,1);
            % V = current(i,j,2);
            
            Ul = current(i,j,1);
            Vb = current(i,j,2);
            Ur = current(i,j,1);
            Vt = current(i,j,2);
            
            dir_swim_speed_single = zeros(4,1);
            
            % Advection velocity, regardless of food
            
            % How fast to swim in each direction, based on percent more food

            per_more_food = reshape(percent_more_food(i,j,:), 4,1);

            % fourWayFishSpeed = mapSpeed( 'linear', per_more_food, fish_speed);
            unit_vector = zeros(size(per_more_food));
            unit_vector( per_more_food > 0 ) = 1;
            fourWayFishSpeed = fish_speed * unit_vector;
            
            % Total speed, current PLUS swimming speed
            % dir_swim_speed = vel + fourWayFishSpeed;
            
            % Fish swim speed versus the current ( net gain or loss )
            vals = [  fish_speed - Vt ...
                      fish_speed + Vb ...
                      fish_speed - Ur ...
                      fish_speed + Ul ];

            
            if happy
                % Swim speeds only where current is stronger than speed
                dir_swim_speed_single( vals < 0 ) = abs ( vals(vals<0) );
                
            else
                % dir_swim_speed_single = [ fourWayFishSpeed(1) + Vt ...
                %                           fourWayFishSpeed(2) + Vb ...
                %                           fourWayFishSpeed(3) + Ur ...
                %                           fourWayFishSpeed(4) + Ul ]
            
                % How fast are they swimming where they WANT to go
                dir_swim_speed_single = [ max( fourWayFishSpeed(1) + Vt, 0 ) ...
                                          max( fourWayFishSpeed(2) - Vb, 0 ) ...
                                          max( fourWayFishSpeed(3) + Ur, 0 ) ...
                                          max( fourWayFishSpeed(4) - Ul, 0 ) ];
                % Add the values where they are being pushed around by the current
                dir_swim_speed_single( vals <= 0 ) = abs ( vals(vals<=0) );
                
            end

            dir_swim_speed_full(i,j,:) = dir_swim_speed_single;

            % dir_swim_speed_single = [ max( fourWayFishSpeed(1) + V, 0 ) ...
            %                    max( fourWayFishSpeed(2) - V, 0 ) ...
            %                    max( fourWayFishSpeed(3) + U, 0 ) ...
            %                    max( fourWayFishSpeed(4) - U, 0 ) ];
            % 
            % dir_swim_speed_full(i,j,:) = dir_swim_speed_single;

        end
    end
end