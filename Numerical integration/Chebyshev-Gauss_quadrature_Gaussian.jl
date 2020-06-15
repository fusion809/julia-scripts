# This approximates
# \int_0^\infty e^(-t^2) dt

# Import package manager
using Pkg;

# Installing required Julia modules
Pkg.add("LinearAlgebra")
using LinearAlgebra;

# Required to use the Airy function and its derivative
Pkg.add("SpecialFunctions")
using SpecialFunctions;

# Number of steps, excluding starting point
N                       = 1000000;
# Our truncated integration domain is [a,b]
a                       = 0.0;
b                       = 1000.0;
# Ai(y) is what we're finding
y                       = 1.0;
# Column vector of integers from 1 to N
n                       = 1:1:N;
# Chebyshev-Gauss grid
tt                      = pi*((2*n.-1)/(2*N));
x                       = cos.(tt);
# Our transformed domain variable
t                       = (b-a)/2.0*x.+(a+b)/2.0;
# Our integrand; sin(tt) is a more efficient way to compute sqrt(1-x^2)
integrand               = sin.(tt).*exp.(-t.^2);
# Chebyshev-Gauss quadrature approximation
Chebyshev_approx        = ((b-a)*pi)/(2*N)*sum(integrand);
# Analytical solution to this problem
exact                   = sqrt(pi)/2;
# Error in our quadrature approximation of the integral
error                   = abs(Chebyshev_approx-exact);

println("error is $(error)")
# Obtains very low errors as our curve is the bell curve and it converges very quickly and does not oscillate. 