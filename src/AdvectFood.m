function [ Food_update ] = AdvectFood( Food, currents, dt, dx, dy)
%ADVECTFOOD Advect the food field with the ocean current for one time step.
%
%   Food_update = AdvectFood(Food, currents, dt, dx, dy)
%
%   Uses first-order upwind differences. The upwind side is chosen once
%   for the whole interior: backward differences in x (y) if every
%   interior u (v) is positive, forward differences otherwise.
%
%   Boundaries: when there is a current in x, the boundary columns are
%   copied from the interior. The first and last rows are copied from
%   their neighbours (zero gradient in y).
%
% INPUTS:
%   Food        : [ny x nx] food concentration [g/m^3]
%   currents    : [ny x nx x 2] current, (:,:,1) = u in x, (:,:,2) = v in y
%                 [m/day]
%   dt          : time step [days]
%   dx, dy      : grid spacing in x and y [m]
%
% OUTPUT:
%   Food_update : [ny x nx] food concentration after one step [g/m^3]
%
% Brzenski, Jared
% 10/16/2026
% -------------------------------------------------------------------------
U = currents(:,:,1); 
V = currents(:,:,2);

% Index maxs
[ n, m ] = size(Food);
% Interior points to iterate over
j = 2:n-1;
i = 2:m-1;
% X(U) derivative holder
dUFood = zeros( n, m );
dUFoodL = zeros( n, m );
dUFoodR = zeros( n, m );
% Y(V) derivative holder
dVFood = zeros( n, m );
dVFoodL = zeros( n, m );
dVFoodR = zeros( n, m );

dUFoodL(j,i) = ( Food(j,i) - Food(j,i-1) ) ./ dx;
dVFoodB(j,i) = ( Food(j,i) - Food(j-1,i) ) ./ dy;
dUFoodR(j,i) = ( Food(j,i+1) - Food(j,i) ) ./ dx;
dVFoodF(j,i) = ( Food(j+1,i) - Food(j,i) ) ./ dy;

% Sided advection based on the sign of apparent velocity
if (U(j,i) > 0)
    dUFood(j,i) = dUFoodL(j,i);
else
    dUFood(j,i) = dUFoodR(j,i);
end

if (V(j,i) > 0)
    dVFood(j,i) = dVFoodB(j,i);
else
    dVFood(j,i) = dVFoodF(j,i);
end

% Update using apparent speed, and derivatives
Food_update = Food - dt.*( U.*dUFood + V.*dVFood);

%% Boundary Conditions
% Periodic everywhere!
% Fish_update = makeFishPositive( Fish_update );
% 
% Fish_update(:,1) = Fish_update(:,end-1); % left
% Fish_update(:,end) = Fish_update(:,2);   % right
% Fish_update(1,:) = Fish_update(end-1,:); % top
% Fish_update(end,:) = Fish_update(2,:);   % bottom
% 
% Fish_update(1,1) = 0;
% Fish_update(end,1) = 0;
% Fish_update(1,end) = 0;
% Fish_update(end,end) = 0;
%Food_update(:,1) = Food_update(:,end-1); % left
%Food_update(:,end) = Food_update(:,2);   % right
%Food_update(1,:) = Food_update(end-1,:); % top
%Food_update(end,:) = Food_update(2,:);   % bottom

% Periodic in X
if max(U(:)) > 0.0
    Food_update(:,end) = Food_update(:,end-1);
    Food_update(:,2) = Food_update(:,end-1);
    Food_update(:,1) = Food_update(:,2);
end
if max(U(:)) < 0.0
    Food_update(:,1) = Food_update(:,2);
    Food_update(:,end) = Food_update(:,end-1);
    Food_update(:,end-1) = Food_update(:,2);
end

% if mean(U(:)) < 1e-8
%     Food_update(:,1) = 0;
%     Food_update(:,end) = 0;
% end

% % Periodic in Y
% if max(V(:)) > 0.0
%     Food_update(1,:) = Food_update(2,:);
%     Food_update(end-1,:) = Food_update(2,:);
%     Food_update(end,:) = Food_update(end-1,:);
% 
% end
% if max(V(:)) < 0.0
%     Food_update(1,:) = Food_update(2,:);
%     Food_update(end-1,:) = Food_update(2,:);
%     Food_update(end,:) = Food_update(end-1,:);
% end

if mean(V(:)) < 1e-8
    Food_update(1,:) = 0;
    Food_update(end,:) = 0;
end

% Neumann in Y
Food_update(1,:) = Food_update(2,:);
Food_update(end,:) = Food_update(end-1,:);

% Check to make sure fish and food is conserved.
%Food_update = conserveQuantity( Food, Food_update);

end
