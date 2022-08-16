function [output, mesh, V] = mom2d_md(mesh, V, N);
%   SYNTAX
%   [output, mesh, V] = mom2d_md(mesh, V, N);
%   DESCRIPTION
%   This function fills MoM matrix and solves MoM equations for an arbitrary 
%   2D transmission-line geometry. It uses an approach given in the paper
%   M. B. Bazdar, A. R. Djordjevic, R. G. Harrington, and T. K. Sarkar, 
%   “Evaluation of quasi-static matrix parameters for multiconductor 
%   transmission lines using Galerkin’s method,” IEEE Trans. Microwave 
%   Theory Tech., vol. 47, no. 7, July 1994, pp. 1223-1228. 
%   for a 2D TL geomtery. Only thick metal conductors are considered. 
%   The charge conservation law is imposed for the total charge. 
%
%   It also employes a simple adaptive mesh refinement routine desribed in
%   particular in A. Bondeson, T.  Rylander, and  P.  Ingelström,  Computational
%   Electromagnetics, Springer, New York, 2005, Series: Texts in Applied Mathematics, 
%   Vol. 51, pp. 153-164, pp. 165-170. The half edge subdivision is used to refine 
%   the mesh in the areas of maximum charge density. 25% mesh refinement means
%   that 25% of edges are subdivided, with the maximum total charge. 
%
%   Inputs:
%   mesh, array of edge parameters including parameters of adjacent dielectric 
%   V, array of voltages on metal conductors 
%   N, number of mesh refinement steps, 25% of edges are divided at every step      
%   N = 1 corresponds to the original mesh
%
%   Outputs:
%   MATLAB structure output with four fields:
%   charge  - total line charge per edge
%   rcond   - conditioning number of the MoM matrix
%   error   - relative error in charge conservation law
%   delta_s - relative error percentage for every refinement step
%   MATLAB structure mesh with refined edges
%   Voltage vector V on refined edges
%   
%   Author: S.N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

const.epsilon       = 8.85418782e-012;      
Mesh                = mesh;     % copy
V_                  = V;        % copy
for nn = 1:N    %   Mesh refinement loop     
    mesh            = Mesh;
    V               = V_; 
    output.charge   = [];
    %------------------------------------------------------------
    % MoM matrix filling (first metal rows, then dielectric rows) 
    nmetal      = find(mesh.e1 == 0);       % Metal edges
    ndiel       = find(mesh.e1 ~= 0);       % Dielectric edges
    Ntotal      = length(nmetal) + length(ndiel);
    Z           = zeros(Ntotal, Ntotal);        % MoM impedance matrix   
    Distance    = zeros(1, Ntotal);             % Distance vector
    Kernel      = zeros(1, Ntotal);             % Integral kernel
    eps         = 1e-6*min(mesh.eh);            % Dummy variable

    % MoM matrix with charge conservation law: divide every metal row by h(m), 
    % subtract the last metal row from the others; and replace the last metal 
    % equation by the total charge conservation law

    % Impedance matrix filling (first nmetal rows)  
    for m = nmetal
        Distance    = sqrt((mesh.ec(1, :) - mesh.ec(1, m)).^2 ...
                         + (mesh.ec(2, :) - mesh.ec(2, m)).^2) + eps;
        Kernel      = log(Distance);
        Z(m, :)     = -1/(2*pi*const.epsilon)*Kernel.*mesh.eh;
        % Diagonal (singular elements)
        Z(m,m)      = -1/(2*pi*const.epsilon)...
                      *mesh.eh(m)*(log(mesh.eh(m))-1.5);    
    end
    for m = nmetal(1:end-1)
        Z(m, :) = Z(m, :) - Z(nmetal(end), :);
        V(:, m) = V(:, m) - V(:, nmetal(end));
    end
    % Last metal row is replaced by the total charge conservation law
    Z(nmetal(end), :)   = ones(1, Ntotal).*mesh.eh/const.epsilon;
    V(:, nmetal(end))   = 0;
    % Impedance matrix filling (last nmetal+1:Ntotal dielectric rows)
    for m = ndiel
        dist_x      = mesh.ec(1, m) - mesh.ec(1, :);
        dist_y      = mesh.ec(2, m) - mesh.ec(2, :);
        dot_pr      = dist_x*mesh.en(1, m) + dist_y.*mesh.en(2, m);    
        Distance    = dist_x.^2 + dist_y.^2 + eps;
        Kernel      = dot_pr./Distance;
        Z(m, :)     = 1/(2*pi*const.epsilon)*Kernel.*mesh.eh*mesh.eh(m)...
                      *(mesh.e1(m) - mesh.e2(m));
    end
    % Diagonal elements
    for m = ndiel
        Z(m,m)      =  -1/(2*const.epsilon)*mesh.eh(m)...
                      .*(mesh.e1(m) + mesh.e2(m));
    end
    % Solution of linear equations (for one or several right-hand sides)
    a       = (Z\V')';      % MoM coefficients a(n) from matrix equation
                            % Batch solution is used
    output.rcond = rcond(Z);
    for m  = 1:size(V, 1)   % Loop over voltage vectors    
        output.charge(m, :) = a(m, :).*mesh.eh;  % Total line charge per edge (C/m)
        output.error(m)     = 100*sum(output.charge(m, :))/max(abs(output.charge(m, :)));
    end
    %   Simple adaptive mesh refinement (error control by static energy)
    %-----------------------------------------------------------
    % find energy (or power) solution error
    for m  = 1:size(V, 1)   % Loop over voltage vectors    
        WV(m)    = 0.5*sum(output.charge.*V(m, :));
    end
    W(nn) = mean(WV);   
    output.delta_s(1) = NaN;
    % stop refinement if N = 1    
    if (N == 1) break; end;       
    if (nn >1)        
        output.delta_s(nn)  = 100*abs(W(nn) - W(nn-1))/abs(W(nn)); % percentage
        semilogy(output.delta_s,  '-.or', 'LineWidth', 2); 
        title('Relative convergence (error percentage)');
        xlabel('iteration number'); ylabel('relative error'); 
        grid on; drawnow; pause(0.25);
    end
    %-----------------------------------------------------------
    % Mesh refinement (so far for one voltage vector only)     
    charge          = abs(output.charge);
    [dummy, index]  = sort(charge);
    index           = index(round(0.75*end):round(1.00*end)); 
                    % 25 percent refinement per step
    % create new mesh - remove old (larger) edges
    Mesh.e(:, index) = [];
    Mesh.ei(index)   = [];
    Mesh.ec(:, index)= [];
    Mesh.eh(index)   = []; 
    Mesh.en(:, index)= []; 
    Mesh.e1(index)   = [];            
    Mesh.e2(index)   = [];
    V_(index)        = [];
    % create new mesh - add new (smaller) edges
    for n = index
        middle_point    = mesh.ec(:, n);
        Mesh.P          = [Mesh.P  middle_point];  % add node in the middle                                          
        node_number     = length(Mesh.P);
        Mesh.e          = [[mesh.e(1, n)  node_number]' ...
                           [mesh.e(2, n)  node_number]' Mesh.e];
        Mesh.ei         = [mesh.ei(1, n) mesh.ei(1, n) Mesh.ei]; 
                                        % add two new ei (the same)
        edge_center1    = 0.5*(mesh.P(:, mesh.e(1, n)) + middle_point);
        edge_center2    = 0.5*(mesh.P(:, mesh.e(2, n)) + middle_point);    
        Mesh.ec         = [edge_center1 edge_center2 Mesh.ec];
        h               = mesh.eh(1, n);                
        Mesh.eh         = [0.5*h 0.5*h Mesh.eh];        
        normal          = mesh.en(:, n);                    
        Mesh.en         = [normal normal Mesh.en ];  % add two new en
        e1              = mesh.e1(1, n);    
        Mesh.e1         = [e1 e1 Mesh.e1];  % add two new e1
        e2              = mesh.e2(1, n);
        Mesh.e2         = [e2 e2 Mesh.e2];
        V_              = [V(1, n) V(1, n) V_]; % add two new V
    end
end
close gcf;
