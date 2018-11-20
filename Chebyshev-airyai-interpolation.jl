# Install PyPlot and SpecialFunctions, if not already installed
using Pkg;
Pkg.add("PyPlot")
Pkg.add("SpecialFunctions")
# Use it
using PyPlot;
using SpecialFunctions;
# Constants
N    = 1000;
NN   = 100000;
OR   = 10;
p0   = -100;
p1   = 10;
dx   = 2/NN;
# Function
function f(x)
    # AiryAi(x) is the function to interpolate
    airyai.(x)
end
# Extrema grid
n    = 0:1:N;
x    = -cos.(pi*n/N);
#y   = cos.(acos.(x)*OR);
# Chebyshev extrema mapped to [p0,p1]
px   = ((p1-p0)/2)*x.+(p1+p0)/2;
y    = f.(px);
# Chebyshev matrix, T[i,j] is Chebyshev T_j(x_i)
T    = cos.(acos.(x)*n');
# Chebyshev coefficients
a    = T\y;
# Linearly spaced grid on [-1,1] with spacing of dx
xx   = -1:dx:1;
# xx mapped to [p0, p1]
pxx  = ((p1-p0)/2)*xx.+(p1+p0)/2;
# List of integers from and including 0 to NN
nn   = 0:1:NN;
# Chebyshev matrix on xx
TT   = cos.(acos.(xx)*n');
# Interpolate f(x) on pxx grid
yy   = TT*a ;
# Error
diff = abs.(airyai.(pxx)-yy);
err  = sqrt(diff'*diff);
# Plot interpolated y
PyPlot.figure(1)
PyPlot.plot(pxx,yy)
# Plot log of error
PyPlot.figure(2)
PyPlot.semilogy(pxx,diff)
