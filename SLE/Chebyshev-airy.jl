# My first Julia program
# Import package manager
using Pkg;

# Install and import ODE
Pkg.add("LinearAlgebra")
using LinearAlgebra;

N         = 1000;
a         = 0;
b         = 10;
k         = 1;
# Column vector of integers from 0 to N
n         = 0:1:N;
nn        = 1:1:N-1;
# Chebyshev extrema grid
x         = -cos.(pi*n/N);
xsub      = x[2:N];
ysub      = (b-a)/2*xsub.+(a+b)/2;
T         = cos.(acos.(x)*n');
Tsub      = T[2:N,:];
sqx       = -((xsub.^2).-1);
Usub      = Diagonal(((sqx).^(-0.5)))*sin.(acos.(xsub)*n');
dTsub     = Usub * Diagonal(n);
dT        = [(-((-1).^n).*n.^2)'; dTsub; (n.^2)'];
D1        = T\dT;
# Second derivative of T(x) on extrema grid without endpoints
d2Tsubr   = Diagonal(xsub) * Usub - Tsub*Diagonal(n);
d2Tsub    = - Diagonal(((xsub.^2).-1).^(-1))*d2Tsubr*Diagonal(n);
# Second derivative of T(x) on extrema grid with endpoints
d2Top     = (((-1).^n).*(n.^2).*((n.^2).-1)/3)';
d2Tbot    = ((n.^2).*((n.^2).-1)/3)';
d2T       = [d2Top; d2Tsub; d2Tbot];
# Second-order differentiation matrix for extrema grid
D2        = D1*D1;
# Second-order differentiation matrix for extrema grid without endpoints
E1        = D1[2:N,2:N];
E2        = D2[2:N,2:N];
H         = -4/((b-a)^2) * E2 + k * Diagonal(ysub);
EIG       = eigen(H);
Y         = EIG.vectors;
Lam       = EIG.values;