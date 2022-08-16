function  D = doping_resistivity(res, flag)
%   SYNTAX
%   D = doping_resistivity(res, flag)
%   DESCRIPTION
%   This function computes doping concentrations (donors or acceptors) based on Si wafer resistivity
%   Inputs:
%   res - resistivity in ohm cm
%   flag - n for donor doping; p for acceptor doping
%   Output:
%   D - doping concentration in cm-3
%   Source:
%   Properties of Crystalline Silicon, Robert Hull, ed., Emis Series, 
%   INSPEC, 1999, pp. 413-416
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

nres(1) =  1.60e-4; pres(1) =  1.30e-4; ND(1) =  1e+21;
nres(2) =  7.70e-4; pres(2) =  1.17e-3; ND(2) =  1e+20;
nres(3) =  5.78e-3; pres(3) =  8.87e-3; ND(3) =  1e+19;
nres(4) =  2.36e-2; pres(4) =  4.35e-2; ND(4) =  1e+18;
nres(5) =  8.38e-2; pres(5) =  2.02e-1; ND(5) =  1e+17;
nres(6) =  5.23e-1; pres(6) =  1.44e+0; ND(6) =  1e+16;
nres(7) =  4.48e+0; pres(7) =  1.33e+1; ND(7) =  1e+15;
nres(8) =  4.29e+1; pres(8) =  1.31e+2; ND(8) =  1e+14;
nres(9) =  4.30e+2; pres(9) =  1.30e+3; ND(9) =  1e+13;
nres(10)=  4.30e+3; pres(10)=  1.30e+4; ND(10)=  1e+12;

if flag=='n'
    D = interp1(log10(nres), log10(ND), log10(res), 'spline');
    D = 10.^D;
end 
if flag=='p'
    D = interp1(log10(pres), log10(ND), log10(res), 'spline');
    D = 10.^D;
end

