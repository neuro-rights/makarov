function [Z, eps_eff] = imp_stripline(w, d, eps_r)
%   SYNTAX
%   [Z, eps_eff] = imp_stripline(w, d, eps_r)
%   DESCRIPTION
%   This function computes characteristic impedance Z and effective diel.
%   constant of a symmetric stripline of zero thickness. The analytical
%   model follows (exact result!) B. Bhat and S. K. Koul, Stripline-Like
%   Transmission Lines for Microwave Integrated Circuits, Wiley, New York,
%   1989, pp.65-66.
%   EXAMPLE: 
%   Tested with
%   http://www.rogerscorporation.com/mwu/mwi_java/Mwij_vp.html
%   w       = 1e-3;
%   d       = 5e-3; NOTE: d is the half -separation distance!
%   eps_r   = 10;
%   [Z eps_eff] = imp_stripline(w, d, eps_r);
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

k1          =   1./cosh(pi*w./(4*d));
k2          =   sqrt(1-k1.*k1);
K1          =   ellipke(k1.*k1);
K2          =   ellipke(k2.*k2);
Z           =   1./(4*const.c*const.epsilon*sqrt(eps_r)).*K1./K2;
eps_eff     = eps_r;


