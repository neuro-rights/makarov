function  [ns, ps, phis, Vbi] = fun_scnm1(ni, NE, VT, x)
%   SYNTAX
%   [ns, ps, phis, Vbi] = fun_scnm1(ni, NE, VT, x)
%   DESCRIPTION
%   This function implements the space-charge neutrality model (SCNM) at
%   zero bias voltage
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

%   the two-domain approach is required for better accuracy:       
ns1   = +0.5*(4*ni^2)./(sqrt(NE.^2 + 4*ni^2) - NE);   
ns2   = +0.5*          (sqrt(NE.^2 + 4*ni^2) + NE);
ns    = [ns1(find(x<=0)) ns2(find(x>0))];   
ps    = ns - NE; 
phis  = VT*log(ns/ni);                          %   when setting Vn = 0    
Vbi   = -(phis(1)-phis(end));                   %   built-in voltage (for test)   