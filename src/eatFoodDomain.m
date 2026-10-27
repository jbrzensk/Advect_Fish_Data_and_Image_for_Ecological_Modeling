function [food_updated] = eatFoodDomain(food, fish, consumption_ratio)
%EATFOODDOMAIN Remove the food eaten by the fish in one time step.
%
%   food_updated = eatFoodDomain(food, fish, consumption_ratio)
%
%   Each cell loses fish .* consumption_ratio .* food, and the result is
%   floored at zero.
%
% INPUTS:
%   food              : [ny x nx] food concentration [g/m^3]
%   fish              : [ny x nx] fish concentration [g/m^3]
%   consumption_ratio : scalar or [ny x nx], fraction of the food eaten per
%                       unit of fish in one step (eating_factor*dt in the
%                       test scripts)
%
% OUTPUT:
%   food_updated      : [ny x nx] food concentration after eating [g/m^3]
%
% NOTES:
%   - Every test script sets eating_factor = 0, so no food is eaten.
%
% Brzenski, Jared
% 02/26/2025
% -------------------------------------------------------------------------
eaten = fish .* consumption_ratio .* food;
food_updated = max(food - eaten, 0);
%
% make sure no negative food
%food_updated( food_updated < 0 ) = 0;
%
end