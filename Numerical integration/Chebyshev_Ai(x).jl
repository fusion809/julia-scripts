# For a reason I cannot figure out, this script gets a
# vastly inaccurate result
# Import package manager
using Pkg;

# Installing required Julia modules
Pkg.add("LinearAlgebra")
using LinearAlgebra;

# Required to use the Airy function and its derivative
Pkg.add("SpecialFunctions")
using SpecialFunctions;

# Number of steps, excluding starting point
N                       = 10000;
# Our truncated integration domain is [a,b]
a                       = 0.0;
b                       = 1000.0;
# Ai(k) is what we're finding
k                       = 1.0;
# Column vector of integers from 0 to N
n                       = 0:1:N;
# Chebyshev extrema grid
t                       = pi*(-n/N.+1.0);
x                       = cos.(t);
# Our transformed domain variable
y                       = (b-a)/2.0*x.+(a+b)/2.0;
# Clearing x as it's unused henceforth
x                       = nothing;
# T_{mn} = T_n(x_m), where T_n is the Chebyshev polynomial of
# the first kind, and x_m are points on the Chebyshev extrema grid
T                       = cos.(t*n');
# T'_n(x_m) = nU_{n-1}(x_m), m=1,2,3,...,N-1
dTsub                   = Diagonal(csc.(t[2:N]))*sin.(t[2:N]*n')*Diagonal(n);
# Unset t
t                       = nothing;
# Adding endpoints, m=0 & N.
dT                      = [(-((-1.0).^n).*n.^2)'; dTsub; (n.^2)'];
dTsub                   = nothing;
# Calculate first order differentiation matrix
D1                      = dT/T;
D1                      = 2/(b-a)*D1;
dT                      = nothing;
T                       = nothing;
integrand               = (cos.(1.0/3.0*y.^3+k*y))/pi;
Ai_Chebyshev_integral   = D1\integrand;
#integrand               = nothing;
Ai_Chebyshev_approx     = Ai_Chebyshev_integral[N+1]-Ai_Chebyshev_integral[1];
Ai_exact                = airyai(k);
error                   = abs(Ai_Chebyshev_approx-Ai_exact);