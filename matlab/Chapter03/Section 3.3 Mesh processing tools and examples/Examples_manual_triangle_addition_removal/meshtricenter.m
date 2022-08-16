function [C] = meshtricenter(P, t) 
%   Triangle centroids (Nx3 array)
%   SNM Winter 2015
    C   = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) + P(t(:, 3), :));
end
   
    
