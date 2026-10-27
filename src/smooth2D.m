function [smoothed] = smooth2D(input,radius)
%SMOOTH2D Mass-preserving box smoothing of a 2D field.
%
%   smoothed = smooth2D(input, radius)
%
%   Replaces each interior point by the mean of the (2*radius+1) x
%   (2*radius+1) box around it. The two outermost rows and columns are
%   replaced by simple neighbour averages. The result is then rescaled so
%   its sum equals the sum of the input.
%
% INPUTS:
%   input    : [ny x nx] field, with ny > 1 and nx > 1
%   radius   : box half-width in grid cells (integer)
%
% OUTPUT:
%   smoothed : [ny x nx] smoothed field with the same sum as input
%
% NOTES:
%   - The test scripts can apply this as artificial diffusion
%     (do_smooth = 1), but it is turned off in every test.
%
% Brzenski, Jared
% 12/17/2024
% -------------------------------------------------------------------------
    [n m] = size(input);
    %
    assert(n>1 && m>1, 'Must give matrix of values, not vector');
    %
    original_total = sum(input, 'all');
    %
    smoothed = input;
    %
    rad_squared = 1 ./ ( (2*radius+1)*(2*radius+1) );
    %
    for ii=2+radius:m-1-radius
        for jj=2+radius:n-1-radius
            holder = 0.0;
            for i=ii-radius:ii+radius
                for j=jj-radius:jj+radius
                    holder = holder + input(j,i);
                end
            end
            smoothed(jj,ii) = holder * rad_squared;
        end
    end
    % Do the edges
    % Left
    for ii=2:n-1
        % 1 away
        smoothed(ii,2) = 0.25*( input(ii,1) + input(ii,3) + input(ii-1,2) + input(ii+1,2) );
        % Edge
        smoothed(ii,1) = 0.5 * (input(ii-1,1) + input(ii,1));
    end
    
    % Right
    for ii=2:n-1
        
        smoothed(ii,end-1) = 0.25*( input(ii,end) + input(ii,end-2) + input(ii-1,end-1) + input(ii+1,end-1) );
        % Edge
        smoothed(ii,end) = 0.5 * (input(ii-1,end) + input(ii,end));
    end
    
    %Bottom
    for jj=2:m-1
        smoothed(2,jj) = 0.25 * ( input(1,jj) + input(3,jj) + input(2,jj-1) + input(2,jj+1));
        smoothed(1,jj) = 0.5 * (input(1,jj-1) + input(1,jj));
    end
     
    %Top
    for jj=2:m-1
        smoothed(end-1,jj) = 0.25 * ( input(end,jj) + input(end-2,jj) + input(end-1,jj-1) + input(end-1,jj+1));
        smoothed(end,jj) = 0.5 * (input(end,jj-1) + input(end,jj) );
    end
    %
    smoothed_total = sum(smoothed, 'all');
    %
    scaling_factor = original_total / smoothed_total;
    %
    smoothed = smoothed * scaling_factor;
    



end
