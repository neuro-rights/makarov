function  [ns, ps, phis, Vbi] = fun_scnm2(ni, NE, VT, V, x)
%   SYNTAX
%   [ns, ps, phis, Vbi] = fun_scnm2(ni, NE, VT, V, x)
%   DESCRIPTION
%   This function implements the space-charge neutrality model (SCNM) at
%   arbitrary bias voltages
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

E     = exp(V/VT);
%   the two-domain approach is required for better accuracy:       
ns1   = +0.5*(4*ni^2*E)./(sqrt(NE.^2 + 4*ni^2*E) - NE);   
ns2   = +0.5*            (sqrt(NE.^2 + 4*ni^2*E) + NE);
ns    = [ns1(find(x<=0)) ns2(find(x>0))];   
ps    = ns - NE; 
phis  = VT*log(ns/ni);                          %   when setting phin = 0    
Vbi   = V -(phis(1)-phis(end));                 %   built-in voltage    