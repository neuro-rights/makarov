function [mesh] = mesh(P, t, re_diel, re_diel_r, set_formula);
%   SYNTAX
%   [mesh] = mesh(P, t, re_diel, re_diel_r, set_formula);
%   DESCRIPTION
%   This function forms MATLAB structure mesh for an arbitrary 2D
%   transmission-line geometry. It finds all nontrivial boundary edges for
%   a 2D triangular mesh P, t based on the data about metal/dielectric
%   domains.
%   
%   Inputs:
%   P, t are arrays of triangle parameters from struct2d;
%   re_diel, re_diel_r are dielectric indicator and dielectric constant 
%   values from geometry text file DIMENSIONS
%
%   Outputs:
%   MATLAB structure mesh with the following fields:
%   mesh.P  - global array of nodes
%   mesh.e  - array of all border edges with nontrivial dielectric contrast
%   mesh.ei - number of geometry subdomain to which the edge belongs
%   mesh.ec - edge center points
%   mesh.eh - edge lengths 
%   mesh.en - edge outer normal vectors (from inner medium 1 to medium 2) 
%   mesh.e1 - dielectric constant from inner side 
%   mesh.e2 - dielectric constant from outer side
%   
%   Author: S.N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

% intialize structure
mesh.P  = P;
mesh.e  = [];
mesh.ei = [];
mesh.ec = [];
mesh.eh = [];
mesh.en = [];
mesh.e1 = [];
mesh.e2 = [];
% collect all border edges for every domain S1, S2, etc. separately; 
% also collect edge parameters
max_domain  = max(t(4, :)); 
count       = 0;
for i = 1:length(re_diel)    
    if (~isempty(findstr(set_formula, num2str(i))) & (count+1<=max_domain))
        count = count + 1;
        ti = t(1:3, find(t(4,:)==count));
        [Edgesm, edges_border, Tp, Tm, Vp, Vm, Vp] = edges(P, ti);        
        edges_i(1, 1:length(Vp))   = count;
        edge_vector         = P(:, edges_border(2, :)) - P(:, edges_border(1, :));
        edge_h              = sqrt(dot(edge_vector, edge_vector, 1));
        edge_center         = (P(:, edges_border(2, :)) + P(:, edges_border(1, :)))/2;               
        edge_normal(1,:)    = -edge_vector(2,:);
        edge_normal(2,:)    = +edge_vector(1,:);
        edge_normal(3,:)    = 0;
        test_vector         = edge_center - P(:, Vp);
        dot_product         = dot(edge_normal, test_vector, 1);
        index               = find(dot_product<0);
        edge_normal(:,index)= - edge_normal(:,index);   % outer normal (from 1 to 2)
        temp                = sqrt(dot(edge_normal, edge_normal, 1));
        edge_normal(1,:)    = edge_normal(1,:)./temp;
        edge_normal(2,:)    = edge_normal(2,:)./temp;
        edge_normal(3,:)    = edge_normal(3,:)./temp;
        e1(1, 1:length(Vp)) = 0;                        % metal by default
        e2(1, 1:length(Vp)) = 1;                        % air by default
        if re_diel(i)
            e1(1,:) = str2num(re_diel_r{i});
        end     
        mesh.e  = [mesh.e edges_border];
        mesh.ei = [mesh.ei edges_i];
        mesh.ec = [mesh.ec edge_center];
        mesh.eh = [mesh.eh edge_h];
        mesh.en = [mesh.en edge_normal];
        mesh.e1 = [mesh.e1 e1];
        mesh.e2 = [mesh.e2 e2];
        clear edges_border edges_i edge_center edge_h edge_normal e1 e2
    end    
end
% find dielectric contrast for metal edges first (all edges are sorted)
MetalEdges  = find(mesh.e1 == 0);
DielEdges   = find(mesh.e1 > 0);
for m = MetalEdges
    temp1 = find(mesh.e(1, DielEdges) == mesh.e(1, m));
    temp2 = find(mesh.e(2, DielEdges) == mesh.e(2, m));
    [dummy] = intersect(temp1, temp2); % common edge number
    if ~isempty(dummy)
        mesh.e2(m) = mesh.e1(DielEdges(dummy));
    end
end
% then remove all diel edges in contact with metal (edges are sorted)
remove = [];
for m = DielEdges
    temp1 = find(mesh.e(1, MetalEdges) == mesh.e(1, m));
    temp2 = find(mesh.e(2, MetalEdges) == mesh.e(2, m));
    [dummy] = intersect(temp1, temp2);  % common edge number
    if ~isempty(dummy)
        remove = [remove m]; 
    end
end
mesh.e(:, remove)  = [];
mesh.ei(:, remove) = [];
mesh.ec(:, remove) = [];
mesh.eh(:, remove) = [];
mesh.en(:, remove) = [];
mesh.e1(:, remove) = [];
mesh.e2(:, remove) = [];
% and after that find contrast for remaining nontrivial diel edges 
MetalEdges  = find(mesh.e1 == 0);
DielEdges   = find(mesh.e1 > 0);
for m = DielEdges
    temp    = DielEdges;
    temp    = setdiff(temp, m);
    temp1 = find(mesh.e(1, temp) == mesh.e(1, m));
    temp2 = find(mesh.e(2, temp) == mesh.e(2, m));    
    [dummy] = intersect(temp1, temp2); % common edge number
    if ~isempty(dummy)
        mesh.e2(m) = mesh.e1(temp(dummy)); 
    end
end
% and after that remove double dielectric edges
remove = [];
for m = DielEdges
    temp    = DielEdges;
    temp    = setdiff(temp, m);
    temp1 = find(mesh.e(1, temp) == mesh.e(1, m));
    temp2 = find(mesh.e(2, temp) == mesh.e(2, m));
    [dummy] = intersect(temp1, temp2);  % common edge number
    if ( (~isempty(dummy)) & (DielEdges(dummy)>m) )
        remove = [remove temp(dummy)]; 
    end
end
mesh.e(:, remove)  = [];
mesh.ei(:, remove) = [];
mesh.ec(:, remove) = [];
mesh.eh(:, remove) = [];
mesh.en(:, remove) = [];
mesh.e1(:, remove) = [];
mesh.e2(:, remove) = [];
% and finally remove trivial dielectric edges (air-to-air)
remove = [];
MetalEdges  = find(mesh.e1 == 0);
DielEdges   = find(mesh.e1 > 0);
for m = DielEdges
    if ((mesh.e2(m)==1)&(mesh.e1(m)==1))
        remove = [remove m];
    end
end
mesh.e (:, remove)  = [];
mesh.ei(:, remove) = [];
mesh.ec(:, remove) = [];
mesh.eh(:, remove) = [];
mesh.en(:, remove) = [];
mesh.e1(:, remove) = [];
mesh.e2(:, remove) = [];