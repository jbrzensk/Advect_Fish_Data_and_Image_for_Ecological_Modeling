function [f_new] = laplacian2D(f_old, dx, dy, kappa)
%LAPLACIAN2D One explicit diffusion step with zero-flux boundaries.
%
%   f_new = laplacian2D(f_old, dx, dy)
%   f_new = laplacian2D(f_old, dx, dy, kappa)
%
%   Computes f_new = f_old + div(kappa * grad(f_old)) with second-order
%   finite differences. Ghost cells copy the boundary values, which gives
%   zero normal gradient (Neumann) at the edges. For scalar kappa this is
%   f_new = f_old + kappa * nabla^2(f_old).
%
%   There is no time step in the formula, so pass kappa already multiplied
%   by dt. The AD test scripts pass kappa*dt_s [m^2].
%
% INPUTS:
%   f_old  : [ny x nx] field, indexed (row = y, column = x)
%   dx     : grid spacing in x (columns) [m]
%   dy     : grid spacing in y (rows) [m]
%   kappa  : diffusion coefficient times the time step [m^2]. Optional,
%            default 1.
%              scalar : uniform diffusion (tested)
%              matrix : [ny x nx] spatially varying diffusion (untested)
%
% OUTPUT:
%   f_new  : [ny x nx] field after one diffusion step
%
% Brzenski, Jared
% 03/09/2026
% -------------------------------------------------------------------------

    if nargin < 4
        kappa = 1; % a safe value, small numbers output, dt stable.
    end

    % Pad f with ghost cells equal to boundary values => zero normal gradient (Neumann BC)
    f_pad = [f_old(:,1), f_old, f_old(:,end)];   % pad left and right columns
    f_pad = [f_pad(1,:); f_pad; f_pad(end,:)];    % pad top and bottom rows

    % Interior indices in the padded array
    ii = 2:size(f_old,1)+1;   % row indices
    jj = 2:size(f_old,2)+1;   % col indices

    if isscalar(kappa)
        % Uniform diffusion: kappa * nabla^2(f)
        d2f_dx2 = (f_pad(ii, jj+1) - 2*f_pad(ii,jj) + f_pad(ii,jj-1)) / (dx^2);
        d2f_dy2 = (f_pad(ii+1, jj) - 2*f_pad(ii,jj) + f_pad(ii-1,jj)) / (dy^2);
        
        %fprintf('Max div_x: %g, Max div_y: %g\n', max(abs(d2f_dx2(:))), max(abs(d2f_dy2(:))));

        f_new = f_old + kappa * (d2f_dx2 + d2f_dy2);

    else % untested
        % Spatially varying kappa: conservative form div(kappa * grad(f))
        % Pad kappa with Neumann BCs
        k_pad = [kappa(:,1), kappa, kappa(:,end)];
        k_pad = [k_pad(1,:); k_pad; k_pad(end,:)];

        % Face-centered kappa: average of the two cells sharing each face
        k_right = 0.5 * (k_pad(ii, jj) + k_pad(ii, jj+1));
        k_left  = 0.5 * (k_pad(ii, jj) + k_pad(ii, jj-1));
        k_up    = 0.5 * (k_pad(ii, jj) + k_pad(ii+1, jj));
        k_down  = 0.5 * (k_pad(ii, jj) + k_pad(ii-1, jj));

        % div(kappa * grad(f)) in each direction
        div_x = (k_right .* (f_pad(ii, jj+1) - f_pad(ii,jj  )) - ...
                 k_left  .* (f_pad(ii, jj  ) - f_pad(ii,jj-1))) / (dx^2);
        div_y = (k_up    .* (f_pad(ii+1, jj) - f_pad(ii,  jj)) - ...
                 k_down  .* (f_pad(ii,   jj) - f_pad(ii-1,jj))) / (dy^2);
        f_new = f_old + div_x + div_y;
    end

end
