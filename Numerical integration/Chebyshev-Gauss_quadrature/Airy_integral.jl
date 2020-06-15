# b has to be relatively small, as with higher b the result gets very inaccurate, likely due to how much our integrand oscillates with increasing t values
# Approximates, using Chebyshev-Gauss quadrature, 
# Ai(y) = \int_0^\infty \cos(t^3/3+y*t) dt
# Import package manager
using Pkg;

# Required to use the Airy function and its derivative
Pkg.add("SpecialFunctions")
using SpecialFunctions;

# Number of steps, excluding starting point
N                       = 1e7;
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
n                       = nothing;
# Our transformed domain variable
t                       = (b-a)/2.0*x.+(a+b)/2.0;
x                       = nothing;
# Our integrand f(x)
integrand               = sin.(tt).*cos.((t.^3)/3.0+y*t);
t                       = nothing;
tt                      = nothing;
# Chebyshev-Gauss quadrature approximation
Chebyshev_approx        = ((b-a))/(2*N)*sum(integrand);
integrand               = nothing;
# Exact solution is Ai(x)
exact                   = airyai(y);
# Absolute difference between Chebyshev approximation and the exact solution
error                   = abs(Chebyshev_approx-exact);

println("error is $(error)")