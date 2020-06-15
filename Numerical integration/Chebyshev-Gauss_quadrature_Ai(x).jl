# b has to be relatively small, as with higher b the result gets very inaccurate, likely due to how much our integrand oscillates
# Import package manager
using Pkg;

# Installing required Julia modules
Pkg.add("LinearAlgebra")
using LinearAlgebra;

# Required to use the Airy function and its derivative
Pkg.add("SpecialFunctions")
using SpecialFunctions;

# Number of steps, excluding starting point
N                       = 100000000;
# Our truncated integration domain is [a,b]
a                       = 0.0;
b                       = 100.0;
# Ai(y) is what we're finding
y                       = 1.0;
# Column vector of integers from 1 to N
n                       = 1:1:N;
# Chebyshev-Gauss grid
tt                      = pi*((2*n.-1)/(2*N));
x                       = cos.(tt);
# Our transformed domain variable
t                       = (b-a)/2.0*x.+(a+b)/2.0;
integrand               = sin.(tt).*cos.((t.^3)/3.0+y*t);
Ai_Chebyshev_approx     = ((b-a))/(2*N)*sum(integrand);
Ai_exact                = airyai(y);
error                   = abs(Ai_Chebyshev_approx-Ai_exact);