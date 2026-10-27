function [ Food ] = makeFood( varargin )
%MAKEFOOD Build an initial food or fish field on the model grid.
%
%   Food = makeFood(scenario, Food, X, Y)
%   Food = makeFood(scenario, Food, X, Y, x0, y0, amplitude, concentration)
%
%   Fills a field the size of Food with one of the scenarios below, then
%   sets the first and last rows and columns to zero. Only the size of the
%   input Food is used; its values are overwritten.
%
%   Scenarios:
%     'gradient'  Linear ramp along the longer grid direction:
%                 1 - X/max(X) if there are more columns than rows,
%                 otherwise Y/max(Y).
%     'slope'     25*X/max(X) if there are more columns than rows,
%                 otherwise Y/max(Y).
%     'hump'      Gaussian of height amplitude centred at (x0, y0),
%                 sigma = 250 km, floored at concentration.
%     'humpfish'  Same as 'hump' with sigma = 15 km. Used for the
%                 initial fish patch.
%
% INPUTS:
%   scenario      : 'gradient', 'slope', 'hump' or 'humpfish'
%   Food          : [n x m] array, sets the size of the output
%   X, Y          : [n x m] grid coordinates from meshgrid [km]
%   x0, y0        : centre of the hump [km]          ('hump', 'humpfish')
%   amplitude     : peak value of the hump            ('hump', 'humpfish')
%   concentration : minimum value outside the hump    ('hump', 'humpfish')
%
% OUTPUT:
%   Food          : [n x m] field, zero on the boundary rows and columns
%
% Brzenski, Jared
% 02/18/2025
% -------------------------------------------------------------------------
numArgs = length(varargin);

switch numArgs
    case 4 % gradient in longer direction [0,1]
        scenario = varargin{1};
        Food = varargin{2};
        X = varargin{3};
        Y = varargin{4};

    case 8 % Gaussian Hump at (x0, y0)
        scenario = varargin{1};
        Food = varargin{2};
        X = varargin{3};
        Y = varargin{4};
        x0 = varargin{5};
        y0 = varargin{6};
        amplitude = varargin{7};
        concentration = varargin{8};

    otherwise
        error(' Unsupported number of args for makeFood');
end

[ n m ] = size( Food );

if strcmp( 'gradient', scenario )
    % Creating gradient (meshgrid) in X and Y
    if m > n % more columns, longer X direction
        scale = max(X(:));
        Food = 1 - X ./ scale .* ones(size(Food));
    else % more rows, longer Y direction
        scale = max(Y(:));
        Food = Y ./ scale .* ones(size(Food));
    end
    Food(:,1)   = 0.0;
    Food(:,end) = 0.0;
    Food(1,:)   = 0.0;
    Food(end,:) = 0.0;

elseif strcmp( 'slope', scenario )
    % Creating gradient (meshgrid) in X and Y
    if m > n % more columns, longer X direction
        scale = max(X(:));
        slope = 25;
        xmid = scale/2;
        Food = slope .* (X ./ scale) .* ones(size(Food));
    else % more rows, longer Y direction
        scale = max(Y(:));
        Food = Y ./ scale .* ones(size(Food));
    end
    Food(:,1)   = 0.0;
    Food(:,end) = 0.0;
    Food(1,:)   = 0.0;
    Food(end,:) = 0.0;

elseif strcmp( 'hump', scenario )
    % Gaussian distribution at (x0, y0), with a fixed concentration outside
    % the hump region
    sigma_x = 250; % Standard deviation in x
    sigma_y = 250; % Standard deviation in y
    initial_concentration = concentration; % average density
    % Define the Gaussian function
    Food = amplitude * exp(-((X-x0).^2 / (2*sigma_x^2) + (Y-y0).^2 / (2*sigma_y^2)));
    Food( Food < initial_concentration) = initial_concentration;

elseif strcmp( 'humpfish', scenario )
    % Gaussian distribution at (x0, y0), with a fixed concentration outside
    % the hump region
    sigma_x = 15; % Standard deviation in x
    sigma_y = 15; % Standard deviation in y
    initial_concentration = concentration; % average density
    % Define the Gaussian function
    Food = amplitude * exp(-((X-x0).^2 / (2*sigma_x^2) + (Y-y0).^2 / (2*sigma_y^2)));
    Food( Food < initial_concentration) = initial_concentration;

else
    error('Could not find scenario to use, returning Food');
end
    Food(:,1)   = 0.0;
    Food(:,end) = 0.0;
    Food(1,:)   = 0.0;
    Food(end,:) = 0.0;

end

        
