function random_weights_full = randomSpeed( apparent_current, fish_speed, percent_more_food_full  )
%RANDOMSPEED Random direction weights for the random advection (RA).
%
%   random_weights_full = randomSpeed(apparent_current, fish_speed, percent_more_food_full)
%
%   For every cell, draws one angle for each of the four directions from a
%   von Mises distribution (circ_vmrnd). The mean angle is that direction
%   (Up = pi/2, Down = 3*pi/2, Right = 0, Left = pi), and the
%   concentration is kappa_max * max(percent_more_food, 0), with
%   kappa_max = 5, so directions without more food get a uniform random
%   angle. Each angle is split into weights on the four directions
%   (decompose_angle), and the four weight vectors are summed.
%
% INPUTS:
%   apparent_current       : [ny x nx x 4] swim speeds; only its size is
%                            used
%   fish_speed             : maximum swim speed; not used
%   percent_more_food_full : [ny x nx x 4] fractional increase in food
%                            toward each neighbour (percentMoreFoodFourWay)
%
% OUTPUT:
%   random_weights_full    : [ny x nx x 4] weights for Up, Down, Right,
%                            Left; they sum to 4 in every cell
%
% NOTES:
%   - Random numbers come from MATLAB's global stream, so call rng first
%     to make a run repeatable.
%
% Brzenski, Jared
% 03/17/2026
% -------------------------------------------------------------------------
    % Directions matrix
    kappa_max = 5; 

    von_misses_angle = [ pi/2, 3*pi/2, 0, pi];

    directions = [
                -1, 0;  % Up
                 1, 0;  % Down
                 0, 1;  % Right
                 0, -1];% Left

    [n, m, ~] = size(apparent_current);

    random_weights_full = zeros(n,m,4);

    for i=1:n
        for j=1:m
            
            app_vel = squeeze(apparent_current(i,j,:) );
            theta_noise = zeros(4,1);
            % 
            % full_velocity_vectors = zeros(4,2);
            % 
            % for d = 1:4
            %     full_velocity_vectors(d, 1) = app_vel(d) * directions(d,1);
            %     full_velocity_vectors(d, 2) = app_vel(d) * directions(d,2);
            % end

            % dir_swim_speed_single = zeros(4,1);
            
            per_more_food = reshape(percent_more_food_full(i,j,:), 4,1);

            % Random drift vector 
            % rand_unit_vector = zeros(size(per_more_food));
            % rand_unit_vector( per_more_food <= 0 ) = 1*fish_speed/10;
            % rand_unit_vector( per_more_food > 0 ) = fish_speed;

            for dd=1:4
                %G = ((per_more_food(dd) + 1) / 2) * kappa_max; % Old 2025
                
                G = max(per_more_food(dd), 0) * kappa_max;
                % G = alpha * abs(per_more_food(dd))
                % theta_noise(dd) = vonmisesrnd(von_misses_angle(dd), G, [1,1]);  % Angular noise per grid cell
                theta_noise(dd) = circ_vmrnd( von_misses_angle(dd), G, 1 );
            end

            % weights = zeros(4,1);
            % up_weight    = decompose_angle( theta_noise(1) );
            % down_weight  = decompose_angle( theta_noise(2) ); 
            % right_weight = decompose_angle( theta_noise(3) );
            % left_weight  = decompose_angle( theta_noise(4) );

            weight = zeros(4,4);
            % Up, Down, Right, Left
            for dd=1:4
                weight(dd,:) = decompose_angle(theta_noise(dd'));
            end

            activeDirections = find( per_more_food > 0); % indices
            totalDirections = numel(activeDirections);   % number of indices
            
            total_weight = zeros(1,4);
            
            for ii=1:4 % totalDirections
                % Going in every direction
                total_weight = total_weight + weight(ii, : );
                % Only in direction of a gradient
                %total_weight = total_weight + weight(activeDirections(ii), : );  
            end
            
            
            random_weights_full(i,j,:) = total_weight;

        end
    end
end

