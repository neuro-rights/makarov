function [Z eps_eff] = imp_slotline_sandwich(w, d, eps_r, f)
%   SYNTAX
%   [Z eps_eff] = imp_slotline_sandwich(w, d, eps_r, f)
%   DESCRIPTION
%   This function computes characteristic impedance Z and effective diel.
%   constant of a sandwich slotline (+/-d).  The analytical model follows
%   Jiri Svacina, "Dispersion Characteristics of Multilayered Slotlines—A
%   Simple Approach," IEEE Trans Microwave Theory Techniques, Vol. 47, no.
%   9, Sep. 1999, pp. 1826-1829
%   EXAMPLE: 
%   Tested with ANSYS HFSS for a shorted resonant slotline. 
%   Bad at low frequencies!  
%   w           = 3e-3;
%   d           = 1e-3;
%   eps_r       = 4;
%   f           = 3.87e9;
%   [Z eps_eff] = imp_slotline_sandwich(w, d, eps_r, f)
%   length      = 2.998e8/(4*sqrt(eps_eff)*f) 
%   r1      = 0.0888*25.4e-3;
%   r2       = 0.2850*25.4e-3;
%   eps_r   = 2.1;  % Teflon
%   [Z eps_eff] = imp_coaxial(r1, r2, eps_r);
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

eps_eff             = 1 + (eps_r - 1)*ellipke(ke_)./ellipke(ke2).*ellipke(k02)./ellipke(k0_); % changed only here (2) 
Z                   = 60*pi./sqrt(eps_eff).*ellipke(k02)./ellipke(k0_);




