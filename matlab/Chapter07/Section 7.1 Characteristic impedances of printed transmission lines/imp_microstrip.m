function [Z eps_eff] = imp_microstrip(w, d, eps_r)
%   SYNTAX
%   [Z eps_eff] = imp_microstrip(w, d, eps_r)
%   DESCRIPTION
%   This function computes characteristic impedance Z and effective diel.
%   constant of a microstrip of zero thickness, width w and height d.
%   EXAMPLE: 
%    w       = 10e-3;
%    d       = 10e-3;
%    eps_r   = 1;
%   [Z eps_eff] = imp_microstrip(w, d, eps_r)
%   
%   Author: S.N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

const.epsilon       = 8.85418782e-012;                  
const.mu            = 1.25663706e-006;                  
const.c             = 1/sqrt(const.epsilon*const.mu);
const.eta           = sqrt(const.mu/const.epsilon);

eps_eff = (eps_r+1)/2 + (eps_r-1)/2./(sqrt(1 + 12*d./w));

if (w./d<1)
    Z = 60/sqrt(eps_eff).*log(8*d./w+w./(4*d));
else
    Z =120*pi/sqrt(eps_eff)./(w./d + 1.393 + 0.667*log(w./d+1.444));
end

