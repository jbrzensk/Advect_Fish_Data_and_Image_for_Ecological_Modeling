function weights = decompose_angle(theta)
%DECOMPOSE_ANGLE Split a direction into weights on the four grid directions.
%
%   weights = decompose_angle(theta)
%
%   Projects the unit vector at angle theta onto Up (0,1), Down (0,-1),
%   Right (1,0) and Left (-1,0), sets negative projections to zero and
%   normalizes the weights to sum to one.
%
% INPUTS:
%   theta   : angle, counterclockwise from the +x axis [rad]
%
% OUTPUT:
%   weights : [4 x 1] non-negative weights for Up, Down, Right, Left,
%             summing to one
%
% NOTES:
%   - These directions use standard axes (Up = +y), not the row/column
%     directions of the swimming code. See the warning in the code.
%
% Brzenski, Jared
% 03/17/2026
% -------------------------------------------------------------------------
    % Unit vector from the angle
    v = [cos(theta), sin(theta)];

    % WARNING: These directions are for converting standard vectors
    % into angles, etc. These are NOT the same as the directions used
    % by the fish swimming logic!!!
    % Directions: Up, Down, Right, Left
    directions = [
         0,  1;  % Up
         0, -1;  % Down
         1,  0;  % Right
        -1,  0]; % Left

    % Compute dot product of v with each direction
    weights = directions * v';  % 4x1 vector of projections

    % Set negative contributions to 0 (no opposition flow)
    weights = max(weights, 0);

    % Optionally normalize so weights sum to 1
    if sum(weights) > 0
        weights = weights / sum(weights);
    end
end
