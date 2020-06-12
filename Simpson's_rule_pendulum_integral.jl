# This script essentially finds the time taken for a
# simple pendulum that starts with zero velocity at the 
# positive x-axis to reach the negative x-axis (theta=-pi)
# using Simpson's rule
# t is time and y is theta 
# Import the required NumPy functions
using Pkg;

Pkg.add("SpecialFunctions")
Pkg.add("QuadGK")
Pkg.add("PyPlot")
using SpecialFunctions;
using QuadGK;
using PyPlot;
# Acceleration rate due to gravity in metres per second
# squared.
g=9.8
# Length of the pendulum in metres.
l=1.0

# Define our integrand
function f(y)
    return 1/sqrt(abs(2.0*g/l*sin(y)))
end

# Number of steps
N=1000000000;
# How close to theta=0 we start our integration and how
# close to theta=-pi we end our intgration.
# At theta=0,-pi our function becomes undefined, hence
# why we cannot start at exact theta=0 or theta=-pi.
tol=1/(25.0*N);
# Integration interval, [y0, yend]
y0=-tol;
yend=-pi+tol;
# Step size
h=(yend-y0)/N;
#
y=y0;

# Our integration function
function Simpson(h,y,i,N)
    if i == 0 || i == N
        return h/3*f(y)
    elseif (i % 2) == 0
        return 2*h/3*f(y)
    else
        return 4*h/3*f(y)
    end
end

t = 0;
# The actual integration
for i=1:N
    t = t + Simpson(h,y,i,N);
    y = y + h
end
# Our SciPy approximation to the integral
T=abs(quadgk(y -> f(y), 0, -pi)[1]);

# Difference between our SciPy approximation and our 
# approximation using Simpson's method
error = abs(T + t)