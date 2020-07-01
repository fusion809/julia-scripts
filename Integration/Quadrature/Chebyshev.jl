# This uses Chebyshev quadrature with the roots grid of T_n(x)
# and weighing function 1/sqrt(1-x^2) to approximate the integral:
# \int_a^b f(x) dx
function chebyshev_quadrature(f, N, a, b)
    n = 1:1:N;
    # Roots grid
    x = -cos.(((2*n.-1)*pi)/(2*N));
    # Linear transformation of the original integration interval [-1,1]
    # to [a,b]
    u = ((b-a)/2)*x.+(a+b)/2;
    int = (((b-a)*pi)/(2*N))*sum((sqrt.((-x.^2).+1)).*f(u));
    return int
end