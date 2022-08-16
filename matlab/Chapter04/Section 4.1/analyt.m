function I = analyt(P)
%   SYNTAX 
%   I = analyt(P) 
%   DESCRIPTION 
%   This function returns double potential MoM self-integrals I = IsIs(1/r)
%   on triangles given by Eq. (4.20) and follows the method described in:
%   Eibert TF, Hansen V. On the calculation of potential integrals for
%   linear source distributions on triangular domains. IEEE Trans Antennas
%   Propag 1995;43 (12):1499–1502
%   Input:
%   P - array of three triangle vertices - rows (a 3x3 matrix)
%   Output:
%   I - double self-integral. To obtain the complete integral, this result
%   should be multiplied by A^2 (triangle area squared)
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

p1=P(1, :);  % three coordinates of the first triangle vertex
p2=P(2, :);  % three coordinates of the second triangle vertex
p3=P(3, :);  % three coordinates of the third triangle vertex

r12=p2-p1;
r23=p3-p2;
r13=p3-p1;

a=sum(r13.*r13);
b=sum(r13.*r23);
c=sum(r23.*r23);
d=a-2*b+c;

A=sqrt(a);
B=sqrt(b);
C=sqrt(c);
D=sqrt(d);

%   Analytical formula for the first integral (Eibert, Hansen, 1995)
%   I(1/r)
N1=(a-b+A*D)*(b+A*C);
D1=(-a+b+A*D)*(-b+A*C);
    
N2=(-b+c+C*D)*(b+A*C);
D2=(b-c+C*D)*(-b+A*C);
    
N3=(a-b+A*D)*(-b+c+C*D);
D3=(b-c+C*D)*(-a+b+A*D);
    
Int     =1/6*(1/A*log(N1/D1) +1/C*log(N2/D2) +1/D*log(N3/D3));  
I       =4 *Int;%  to obtain the complete integral, 
                %  this result should be multiplied by A^2 (triangle area squared)