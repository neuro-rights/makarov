function [Z, eps_eff,eps_efft] = imp_cpw(w, d, S, t, eps_r);
%   SYNTAX
%   [Z, eps_eff,eps_efft] = imp_cpw(w, d, S, t, eps_r);
%   DESCRIPTION
%   This function computes characteristic impedance, effective dielectric
%   constant, and the conductor-thickness dependent effective dielectric
%   constant for a coplanar waveguide. The analytical model follows
%   "Antenna Engineering Handbook", John L. Volakis, (Mcgraw Hill 4th Ed).
%   The Z-espression is derived from "Transmission line Design Handbook" by
%   B.C.Wadell.
% Inputs :      w - center conductor width
%               d - dielectric substrate thickness
%               S - slot width from adjacent ground planes 
%               t - conductor thickness
%               eps_r - relative dielectric constant
% Outputs:     Z  - characteristic impedance
%              eps_eff - effective dielectric constant
%              eps_efft- conductor-thickness dependent effective dielectric
%   EXAMPLE: 
%   To verify the formulas an online calculator was used at
%   http://chemandy.com/calculators/coplanar-waveguide-with-ground-calculator.htm
%   Note: - 1.) They don't have an option for conductor thickness
%           2.) Website uses S for conductor thickness and w for gap width
%           (notations interchanged).
%   Usage:
%     w       = 3;
%     d       = 4;
%     S       = 1;
%     t       = 0.1;
%     eps_r   = 4.2;
%     [Z, eps_eff,eps_efft] = imp_cpw(w, d, S, t, eps_r);
%     Z_web       = 61.33
%     eps_eff_web = 2.699
%   
%   Authors: V. Iyer, S.N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

% Calculate all relevant parameters
a           = w;
b           = w + 2*S;
a_t         = a + (1.25.*t.*(1 + log(4*pi*a./t))./ pi);
b_t         = b - (1.25.*t.*(1 + log(4*pi*a./t))./ pi);

k           = a./b;
kprime      = sqrt(1 - k.^2);

k_t         = a_t./b_t;
kprime_t    = sqrt(1 - k_t.^2);

k_1         = sinh(pi.*a_t./4./d)./sinh(pi.*b_t./4./d);
kprime_1    = sqrt(1 - k_1.^2);

% Calculate eps_eff
eps_eff     = 1 + ((eps_r -1).*ellipke(kprime).*ellipke(k_1)./2./ellipke(k)./ellipke(kprime_1));
%Calculate eps_efft
p           = (b - a).*ellipke(k)./ 1.4./t./ellipke(kprime);
eps_efft    = eps_eff - ( (eps_eff-1)./ (p+1));

%Calculate Z
Z           = 30.*pi.*ellipke(kprime_t)./sqrt(eps_efft)./ellipke(k_t); 



