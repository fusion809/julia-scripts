# This script essentially finds the time taken for a
# simple pendulum that starts with zero velocity at the 
# positive x-axis to reach the negative x-axis (theta=-pi)
# using Simpson's rule
# t is time and y is theta
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
    # abs is used below to prevent complex number errors due to the sqrt
    return 1/sqrt(abs(2.0*g/l*sin(y)))
end

# Number of steps
N=1000000000;
# How close to theta=0 we start our integration and how
# close to theta=-pi we end our intgration.
# At theta=0,-pi our function becomes undefined, hence
# why we cannot start at exactly theta=0 or end at exactly
# theta=-pi.
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
    if i == 1 || i == N+1
        return h/3*f(y)
    elseif (i % 2) == 1
        return 2*h/3*f(y)
    else
        return 4*h/3*f(y)
    end
end

# Function is used for the integration because otherwise
# we have to use global t within the loop to prevent the 
# error mentioned here (
# https://discourse.julialang.org/t/i-define-t-then-in-a-loop-im-told-its-undefined-why/41258/5)
# and that is bad for performance. Apparently in Julia 1.5, we 
# won't need to do this anymore.

function simpson(y0)
    t = 0.0;
    y = y0;
    # The actual integration
    for i=1:N+1
        t = t + Simpson(h,y,i,N);
        if i < N+1
            y = y + h;
        end
    end

    return t
end

t=simpson(y0)
# Our QuadGK approximation to the integral
T=abs(quadgk(y -> f(y), 0, -pi)[1]);

# Difference between our QuadGK approximation and our 
# approximation using Simpson's method
error = abs(T + t)