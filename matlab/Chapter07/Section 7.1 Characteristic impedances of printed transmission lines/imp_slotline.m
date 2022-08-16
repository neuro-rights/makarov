function [Z eps_eff] = imp_slotline(w, d, eps_r, f)
%   SYNTAX
%   [Z eps_eff] = imp_slotline(w, d, eps_r, f)
%   DESCRIPTION
%   This function computes characteristic impedance Z and effective diel.
%   constant of a simple unshielded slotline (0-d). The analytical model
%   follows Jiri Svacina, "Dispersion Characteristics of Multilayered
%   Slotlines—A Simple Approach," IEEE Trans Microwave Theory Techniques,
%   Vol. 47, no. 9, Sep. 1999, pp. 1826-1829
%   EXAMPLE: 
%   Tested with
%   http://www1.sphere.ne.jp/i-lab/ilab/tool/sl_line_e.htm
%   w       = 2e-3;
%   d       = 0.5e-3;
%   eps_r   = 5;
%   f       = 10e9;
%   [Z eps_eff] = imp_slotline(w, d, eps_r, f)
%   
%   Author: S. N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.


const.epsilon       = 8.85418782e-012;                  
const.mu            = 1.25663706e-006;                  
const.c             = 1/sqrt(const.epsilon*const.mu);
const.eta           = sqrt(const.mu/const.epsilon);

alpha               = tanh(pi*w./(2*d));
ke2                 = 2*alpha./(1 + alpha);
ke_                 = 1 - ke2;

h0                  = d*(1 + 0.01333/(eps_r + 2)*(const.c./f./d).^2); 
alpha0              = tanh(pi*w./(2*h0));
k02                 = 2*alpha0./(1+alpha0);
k0_                 = 1-k02;

eps_eff             = 1 + (eps_r-1)/2.*ellipke(ke_)./ellipke(ke2).*ellipke(k02)./ellipke(k0_);
Z                   = 60*pi./sqrt(eps_eff).*ellipke(k02)./ellipke(k0_);




