function [Z_e, Z_o, eps_eff_e, eps_eff_o] = imp_microstrip_coupled(w, d, S, eps_r);
%   SYNTAX
%   [Z_e, Z_o, eps_eff_e, eps_eff_o] = imp_microstrip_coupled(w, d, S, eps_r);
%   DESCRIPTION
%   This function computes characteristic impedance/effective diel.
%   constant (no dispersion/attenuation) for a coupled microstrip.
%   The analytical model follows
%   E. Hammerstad and O. Jensen, “Accurate models for microstrip
%   computer-aided design,” in IEEE MTT-S Inter. Microwave Symp. Dig.,
%   Washington, DC, 1980, pp. 407-409.
%
%   See also:
%   M. Kirschning and R. H. Jansen, “Accurate wide-range design equations
%   for the frequency-dependent characteristic of parallel coupled
%   microstrip lines,” IEEE Trans. Microwave Theory and Techniques, vol.
%   MTT-32, pp. 83-89, 1984, with corrections in vol. MTT-33, p. 288, 1985.
%   (*) Ch. Wan, “Analytically and accurately determined quasi-static
%   parameters of coupled microstrip lines,” IEEE Trans. Microwave Theory
%   and Techniques, vol. 44, no. 1, Jan. 1996,  pp. 75 – 80.
%   EXAMPLE: 
%   Tested with Ref. (*).   
%   d               = 1.57e-3;       %   Microstrip height, m
%   w               = 0.5*d;         %   Microstrip width, m      
%   S               = 0.2*d;         %   Trace separation, m  
%   eps_r           = 10.0;          %   Relative dielectric constant   
%   [Z_e, Z_o, eps_eff_e, eps_eff_o] = imp_microstrip_coupled(w, d, S, eps_r);
%   Z_e
%   Z_o
%   eps_eff_e
%   eps_eff_o
%   %   Bryant & Weiss (Table I of Ref. (*))
%   data.Z_e_test        = 90.40;
%   data.Z_o_test        = 39.15;
%   data.eps_eff_e_test  = 6.77;
%   data.eps_eff_o_test  = 5.57;
%   data
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

u       = w/d;
g       = S/d;

%   Coefficients
alpha = 0.5*exp(-g);

beta  = 0.2306 + 1/301.8*log(g.^10./(1+(g/3.73).^10)) + 1/5.3*log(1+0.646*g.^1.175);

phi   = 0.8645*u.^0.172;

psi   = 1 + g/1.45 +g.^2.09/3.95;

theta = 1.729+1.175*log(1 + 0.627./(g+0.327*g.^2.17)); 

m     = 0.2175 + (4.113 + (20.36./g).^6 ).^(-0.251) +1/323*log(g.^10./(1+(g./13.8).^10));

n     = log( (10+68.3*g.^2)./(1+32.5*g.^3.093) ).*...
        (1/17.7 + exp( -6.424-0.76*log(g)-(g/0.23).^5) );

mu    = g.*exp(-g) + u.*(20+g.^2)./(10+g.^2);

r     = 1 + 0.15*(1 - exp(1-(eps_r-1).^2/8.2)./(1+g.^-6) );

fo1   = 1 - exp( -0.179*g.^0.15-0.328*g.^r./log(exp(1)+(g/7).^2.8) );

p     = exp(-0.745*g.^0.295)./cosh(g.^0.68);

q     = exp(-1.366-g);

fo    = fo1.*exp( p.*log(u) + q.*sin(pi*log(u)/log(10)) );

phi_e = phi./( psi.*(alpha.*u.^m + (1-alpha).*u.^(-m)) );

phi_o = phi_e -(theta./psi).*exp( beta.*u.^(-n).*log(u) );

% Prop. constants/eff. epsilons

a    = 1 + 1/49*log( (u.^4+(u/52).^2)/(u.^4+0.432) ) ...
         + 1/18.7*log(1+(u/18.1).^3);
     
b    = 0.564*( (eps_r-0.9)./(eps_r+3) ).^0.053;
     
F_o   = fo.*(1+10./u).^(-a.*b); 

F_e   = (1 + 10./mu).^(-a.*b);

eps_eff_e = (eps_r + 1)/2 + (eps_r -1)/2.*F_e;

eps_eff_o = (eps_r + 1)/2 + (eps_r -1)/2.*F_o;

%   Even/odd mode impedances 

fu  = 6 + (2*pi-6)*exp( -(30.666/u).^0.7528 );

Z_0 = const.eta/(2*pi).*log(fu./u + sqrt(1+(2./u).^2));

Z_e = Z_0./(1 - Z_0.*phi_e/const.eta)./sqrt(eps_eff_e);

Z_o = Z_0./(1 - Z_0.*phi_o/const.eta)./sqrt(eps_eff_o);
