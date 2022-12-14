function [A, H] = staticsh(a, b, center, r0, N, nx, ny, nz)
%   SYNTAX
%   [A, H] = staticsh(a, b, center, r0, N, nx, ny, nz)
%   DESCRIPTION
%   STATICSH finds magnetic field generated by an arbitrary ellipse of a
%   uniform current
%   Inputs:
%   a - major axis, m
%   b - minor axis, m
%   loop current of 1A is assumed
%   center - lop center (a 1 by 3 vector)
%   r0 - observation point for the B-field, m (a N by 3 vector) 
%   N  - number of integration points
%   nx, ny, nz - normal vector of the ellipse. Default: (0 0 1)
%   Outputs:
%   A -  magnetic vector potential (N, 3), T*m
%   H -  magnetic field density (N, 3), A/m
%
%   By Sheila Werth, Kaung Myat Win, and S. Makarov 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

    %    EM DATA
    mu0 = 1.25663706e-006;    %    magnetic permeability of vacuum(~air)

    %   NORMALIZE THE NORMAL VECTOR
    temp    = sqrt(nx^2 + ny^2 + nz^2);
    nx      = nx/temp;
    ny      = ny/temp;
    nz      = nz/temp;
    
    %   DISCRETIZATION
    t = linspace(0, 2*pi, N);  %   parameteric ellipse form
    dt = t(2) - t(1); t = t + dt/2; t(end) = [];
    %   INTEGRATION
    H = zeros(size(r0, 1), 3);  %   magnetic field
    A = zeros(size(r0, 1), 3);  %   magnetic vector potential
    rshift0 = r0 - repmat(center, [size(r0, 1) 1]);
    
    a1 = a;    
    b1 = b;   
    
    for m = 1:length(t)
        x  = [+a1*cos(t(m)) +b1*sin(t(m)) 0];       %   parameterization
        dl = [-a1*sin(t(m)) +b1*cos(t(m)) 0]*dt;    %   parameterization
        x  = meshrotate(x,  [-ny nx 0], acos(nz));
        dl = meshrotate(dl, [-ny nx 0], acos(nz));
        dL = repmat(dl, [size(r0, 1) 1]);
        r  = rshift0 - repmat(x, [size(r0, 1) 1]); 
        R  = dot(r, r, 2);
        R1 = R.^0.5;
        R3 = R.^1.5;
        R1 = repmat(R1, [1, 3]);
        R3 = repmat(R3, [1, 3]);        
        H  = H + cross(dL, r, 2)./R3;
        A  = A + dL./R1;
    end
    H = H/(4*pi);
    A = mu0*A/(4*pi);    
end
