function [t] = meshsegtrintersection(orig0, dir0, dist0, P0, t0)
% Code based on the MT algorithm:
% Tomas Moeller and Ben Trumbore, “Fast, Minimum Storage Ray/Triangle Intersection”, 
% Journal of Graphics Tools, 2(1):21—28, 1997
% Authors:
% Janakinadh Yanamadala (jyanamadala@wpi.edu, Vishal Rathi (vkrathi@wpi.edu)
% SNM (makarov@wpi.edu)
% Summer 2014
% Output:
% t -  Distances of point of intersection from the origin of the segment
% (N x 1). If there is no intersection, the corresponding field is zero.  
%(u,v) - Barycentric coordinates for the point of intersection (optional)
% Input:
% orig0 - Origin of the segment (1 x 3)
% dir0 - Normalized direction of the segment from origin (1 x 3)
% dist0 - Length of the segment (1 x 1)
% P0, t0 - triangulation to be tested 

    vert1 = P0(t0(:, 1),:);
    vert2 = P0(t0(:, 2),:);
    vert3 = P0(t0(:, 3),:);
    orig = repmat(orig0, size(vert1, 1),1);             
    dist = repmat(dist0, size(vert1, 1),1);
    dir  = repmat(dir0, size(vert1, 1),1);

    % Initialization of u,v and t
    u = zeros (size(vert1,1),1);
    t = u; v = u;

    % Finding edges
    edge1 = vert2 - vert1;
    edge2 = vert3 - vert1;

    tvec = orig - vert1;                            % distance to vert1 from segment origin
    pvec = cross(dir, edge2, 2);                    % parameter to calculate u
    det  = dot(edge1, pvec, 2);                     % Determinant of matrix M
    parallel = abs(det)< 1024*eps*max(abs(det));    % To test parallel edges with the segment
    if all(parallel)                                % if all parallel than no intersections
        return;
    end

    det(parallel) = 1;              % To avoid division by zero
    inv_det = 1.0 ./ det;           % Finding inverse of the determinant
    u = dot(tvec,pvec,2);           % Calculate u parameter
    u = u.*inv_det;

    % Conditional tests for u and v
    layer1 = (~ parallel & u<0 | u>1);
    if all(layer1)
        return;
    end

    qvec (~layer1,:) = cross(tvec(~layer1,:), edge1(~layer1,:), 2);             % Parameter to calculate v
    v (~layer1,:) = dot(dir(~layer1,:),qvec(~layer1,:),2).*inv_det(~layer1,:);  % Calculate v
    layer2 = (v<=0 | u+v>1);
    if all(layer2)
        return;
    end

    layer = (~layer1&~layer2);
    t(layer,:) = dot(edge2(layer,:),qvec(layer,:),2).*inv_det(layer,:);         % Calculate t
    t(t<0 | t>dist) = 0;                                                        % Comparing distance and t
    t(parallel) = 0;                                                            % Avoiding values of t in parallel cases
    t(isnan(t))= 0;                                                             % Avoiding NaN (Not-a-Number) when right angled triangles are present
end


