function eddy = eddyatpoint(r0, z0, r, z, I0, mu0, omega, sigma)
%   SYNTAX
%   eddy = eddyatpoint(r0, z0, r, z, I0, mu0, omega, sigma)
%   DESCRIPTION
%   Analytical solution for the eddy current density in the conducting
%   material given by Eqs. 11.65 and 11.66 at a point
%    
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

    %   General integration variables
    NMIN = 1000;                                        %   Minimum number of integration points
    STOP = 1e-3;                                        %   Min value of the integrand to be tracked (compared to the max value)
    stop        = 1e-3;                                 %   Min value of the exponential factor exp(alpha(l-z)) to be integrated
    limit       = log(1/stop)/(-z+z0);                  %   Integration limit for the entire integral found from this criterion
    T           = 2*pi/min(r, r0);                      %   Minimum period of oscillating Bessel functions (min 1,000 nodes)
    step1       = min(T/10, limit/NMIN);                %   Step size to resolve Bessel function oscillations
    
    alpha       = [0+step1/2:step1:limit-step1/2];      %   Integration variable alpha 
    N           = length(alpha);                        %   Number of integration points
    beta        = sqrt(alpha.^2 + j*omega*mu0*sigma);   %   Parameter beta
    array1      = besselj(1, alpha*r0);                 %   Integration array
    array2      = besselj(1, alpha*r);                  %   Integration array
    array3      = exp(beta*z-alpha*z0).*...
                                2.*alpha./(alpha+beta); %   Integration array
    %   Integral
    Integral3   = step1*sum(array1.*array2.*array3);
    eddy        = mu0*I0*r0/2*Integral3;
    eddy        = -j*omega*sigma*eddy; 
end
