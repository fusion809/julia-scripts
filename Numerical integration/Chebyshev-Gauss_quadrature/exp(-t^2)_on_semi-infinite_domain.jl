# This approximates
# \int_0^\infty e^{-t^2} dt
# (chosen as it is known it rapidly converges, hence 
# approximation using quadrature should be relatively
# easy).
# Using Chebyshev-Gauss quadrature

# Number of steps, excluding starting point
N                       = 2e8;
# Our truncated integration domain is [a,b]
a                       = 0.0;
b                       = 10.0;
# Column vector of integers from 1 to N
n                       = 1:1:N;
# Chebyshev-Gauss grid
tt                      = pi*((2.0*n.-1.0)/(2.0*N));
x                       = cos.(tt);
# Our transformed domain variable
t                       = (b-a)/2.0*x.+(a+b)/2.0;
# Free up RAM by clearing x
x                       = nothing;
# Our integrand f(x); sin(tt) is a more efficient way to compute sqrt(1-x^2)
integrand               = sin.(tt).*exp.(-t.^2);
# Free up RAM by clearing tt, n and t
tt                      = nothing;
n                       = nothing;
t                       = nothing;
# Chebyshev-Gauss quadrature approximation
Chebyshev_approx        = ((b-a)*pi)/(2.0*N)*sum(integrand);
# Free up RAM by clearing integrand
integrand               = nothing;
# Analytical solution to this problem
exact                   = sqrt(pi)/2.0;
# Error in our quadrature approximation of the integral
error                   = abs(Chebyshev_approx-exact);

println("error is $(error)")
# Obtains very low errors as our curve is the bell curve and it converges very quickly and does not oscillate. 