function [Z eps_eff] = imp_coaxial(r1, r2, eps_r)
%   SYNTAX
%   [Z eps_eff] = imp_coaxial(r1, r2, eps_r);
%   DESCRIPTION
%   This function computes characteristic impedance Z and effective diel.
%   constant of a coaxial line. The analytical model follows
%   D. Pozar, Microwave Engineering, Wiley, New York, 2005, pp. 55-57.
%   EXAMPLE: 
%   (for RG214/U see Pozar, p. 689):
%   r1      = 0.0888*25.4e-3;
%   r2       = 0.2850*25.4e-3;
%   eps_r   = 2.1;  % Teflon
%   [Z eps_eff] = imp_coaxial(r1, r2, eps_r);
%   
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

const.epsilon       = 8.85418782e-012;                  
const.mu            = 1.25663706e-006;                  
const.c             = 1/sqrt(const.epsilon*const.mu);
const.eta           = sqrt(const.mu/const.epsilon);

eps_eff = eps_r;

Z = sqrt(const.mu./(const.epsilon*eps_r)).*log(r2./r1)/(2*pi);
