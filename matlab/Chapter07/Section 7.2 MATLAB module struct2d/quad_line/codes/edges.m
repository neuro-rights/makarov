function [Edgesm, Edgesb, Tp, Tm, Vp, Vm, Vb] = edges(P, t)
%   SYNTAX
%   Edgesm, Edgesb, Tp, Tm, Vp, Vm, Vb] = edges(P, t)
%   DESCRIPTION
%   This function finds the inner and outer edges of a triangular mesh.  For the
%   inner edges, EDGES also finds the two triangles adjacent to that edge
%   and the opposite vertices of each of those triangles.  If an edge is
%   adjacent to more than two triangles, EDGES attempts to avoid returning
%   a (Tp, Tm) pair in the same horizontal plane as each other.
%   
%   Inputs:
%   P, array of nodes
%   t, array of triangles
%
%   Outputs:
%   Edgesm - array of inner edges
%   Edgesb - array of boundary (outer) edges
%   Tp - array of adjacent "plus"  triangles for every inner edge
%   Tm - array of adjacent "minus" triangles for every inner edge
%   Vp - array of opposite vertices for adjacent "plus"  triangles
%   Vm - array of opposite vertices for adjacent "minus" triangles
%   Vb - array of opposite vertexes for all border edges
%   
%   Authors: A. Marut, S.N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

[Edges TpTm]    = findfaces(t(1:3,:), 'edges'); %   Each edge occurs once for each bounded triangle
                                                %   TpTm - list of triangles for each edge
[Edges I]       = sortrows(Edges');             %   All occurrences of each edge are contiguous  
Edges           = Edges';
TpTm            = TpTm(I);

n       = size(Edges,2);
output  = zeros(6,n);
% Rows 1 + 2: Edges, represented by endpoints
% Row 3: Tp, represented by an index into t
% Row 4: Tm, represented by an index into t
% Row 5: Vp, represented by a row index into t (column index is Tp)
% Row 6: Vm, represented by a row index into t (column index is Tn)

i = 1;
triangles_current_edge = [];
vertices_current_edge = [];
vertices_current_edge_indices = [];

for j = 1:n+1 % the purpose of n=1 is to do final processing
    if j <= n
        E  = Edges(:,j);
    end
    
    if j==1 || (j ~= n+1 && all(E == last))
        % make a list of all triangles and opposite vertices for the current edge
        triangles_current_edge(end+1) = TpTm(j);
        vertices_current_edge(end+1)  = sum(t(1:3,triangles_current_edge(end))) - sum(E);
        vertices_current_edge_indices(end+1) = find(t(1:3,triangles_current_edge(end)) == vertices_current_edge(end));
    else
        % process the aforementioned list        
        if length(triangles_current_edge) == 1
            % do nothing
        elseif length(triangles_current_edge) == 2
            output(:,i) = [last; triangles_current_edge'; vertices_current_edge_indices'];
            i = i+1;
        else
            N = [];
            for m = 1:length(vertices_current_edge)
                %Vertexes    = P(:, t(:, m));
                %C1          = Vertexes(:,1) -Vertexes(:,2);
                C1 = P(:,last(1)) - P(:,vertices_current_edge(m));
                C2 = P(:,last(2)) - P(:,vertices_current_edge(m)); % try to vectorize this later
                %C2          = Vertexes(:,1) -Vertexes(:,3);
                
                N(:, m)     = cross(C1, C2);
                N(:, m)     = N(:,m)/norm(N(:,m));    
            end
            
            % assume size(N,2) == 3
            D1 = dot(N(:,2),N(:,3));
            D2 = dot(N(:,3),N(:,1));
            D3 = dot(N(:,1),N(:,2));
            DMin = min([D1, D2, D3]);
            
            % store the results in output
            vertices_current_edge = vertices_current_edge_indices;
            if D1 == DMin
                output(:,i)   = [last; triangles_current_edge(2); triangles_current_edge(1); vertices_current_edge(2); vertices_current_edge(1)];
                output(:,i+1) = [last; triangles_current_edge(3); triangles_current_edge(1); vertices_current_edge(3); vertices_current_edge(1)];
            elseif D2 == DMin
                output(:,i)   = [last; triangles_current_edge(3); triangles_current_edge(2); vertices_current_edge(3); vertices_current_edge(2)];
                output(:,i+1) = [last; triangles_current_edge(1); triangles_current_edge(2); vertices_current_edge(1); vertices_current_edge(2)];
            else
                output(:,i)   = [last; triangles_current_edge(1); triangles_current_edge(3); vertices_current_edge(1); vertices_current_edge(3)];
                output(:,i+1) = [last; triangles_current_edge(2); triangles_current_edge(3); vertices_current_edge(2); vertices_current_edge(3)];
            end
            i = i + 2;
            
        end
        
        % clear the lists and start storing data for the new edge
        if j ~= n+1
            triangles_current_edge = TpTm(j);
            vertices_current_edge  = sum(t(1:3,triangles_current_edge(end))) - sum(E);
            vertices_current_edge_indices = find(t(1:3,triangles_current_edge) == vertices_current_edge);
        end
    end

    last = E;
    
end

output = output(:, 1:i-1);
Edgesm = output(1:2, :);
Edgesb = setdiff(Edges',Edgesm','rows')';
Tp = output(3, :);
Tm = output(4, :);
Vp = output(5, :);
Vm = output(6, :);
Vb = zeros(1, length(Edgesb));

for m = 1:length(Edgesb)
    temp1 = find(t(1, :)==Edgesb(1, m)); % triangle numbers with vertex 1 
    temp2 = find(t(2, :)==Edgesb(2, m)); % triangle numbers with vertex 2 
    [dummy ] = intersect(temp1, temp2); % common triangle number    
    if ~isempty(dummy)
        Vb(m) = t(3, dummy); 
    end
    temp1 = find(t(2, :)==Edgesb(1, m)); % triangle numbers with vertex 1 
    temp2 = find(t(3, :)==Edgesb(2, m)); % triangle numbers with vertex 2 
    [dummy ] = intersect(temp1, temp2); % common triangle number    
    if ~isempty(dummy)
        Vb(m) = t(1, dummy); 
    end
    temp1 = find(t(1, :)==Edgesb(1, m)); % triangle numbers with vertex 1 
    temp2 = find(t(3, :)==Edgesb(2, m)); % triangle numbers with vertex 2 
    [dummy ] = intersect(temp1, temp2); % common triangle number    
    if ~isempty(dummy)
        Vb(m) = t(2, dummy); 
    end    
end

    
    
    
