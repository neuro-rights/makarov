function [Z_e, Z_o] = imp_stripline_coupled(w, d, S, eps_r)
%   SYNTAX
%   [Z_e, Z_o] = imp_stripline_coupled(w, d, S, eps_r)
%   DESCRIPTION
%   This function computes characteristic impedance Z (even and odd modes)
%   of a coupled stripline of zero thickness. The analytical model follows
%   S. B. Cohn, “Shielded coupled strip transmission line,” IRE Trans., 
%   vol. MTT-5, Oct. 1955, pp. 29-37.
%   EXAMPLE: 
%   w       = 2.66e-3;             % Strip width (m)
%   d       = 3.20e-3;             % Half (!) separation distance between ground plates(m)
%   S       = 1e-3;                % Separation between conductors (m)
%   eps_r   = 2.2;                 % Relative dielectric constant
%
%   Authors: Ali Ilhan/S.N. Makarov
%   
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

const.epsilon       = 8.85418782e-012;                  
const.mu            = 1.25663706e-006;                  
const.c             = 1/sqrt(const.epsilon*const.mu);
const.eta           = sqrt(const.mu/const.epsilon); 

b = 2*d;
%for even mode assuming zero conductor thickness
ke          =   tanh(pi/2*w./b).*tanh(pi/2*(w+S)./b);
ke_pr       =   sqrt(1-ke.*ke);
K1          =   ellipke(ke_pr.*ke_pr);
K2          =   ellipke(ke.*ke);
Z_e         =   30*pi./(sqrt(eps_r)).*K1./K2;

% Cf_0 = 2/pi*log(2)*const.epsilon;
% Cf_e = 2/pi*log(1+tanh(pi*S/(2*b)))*const.epsilon;
% Z_e = 94.15/(sqrt(eps_r)*(w/b+1/(2*const.epsilon)*(Cf_0+Cf_e)));

%for odd mode assuming zero conductor thickness
ke          =   tanh(pi/2*w./b)./tanh(pi/2*(w+S)./b);
ke_pr       =   sqrt(1-ke.*ke);
K1          =   ellipke(ke_pr.*ke_pr);
K2          =   ellipke(ke.*ke);
Z_o         =   30*pi./(sqrt(eps_r)).*K1./K2;

% Cf_0 = 2/pi*log(2)*const.epsilon;
% Cf_o = 2/pi*log(1+coth(pi*S/(2*b)))*const.epsilon;
% Z_o = 94.15/(sqrt(eps_r)*(w/b+1/(2*const.epsilon)*(Cf_0+Cf_o)));