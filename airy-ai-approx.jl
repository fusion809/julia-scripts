# This uses the equation:
# Ai(x) = int_0^inf cos(t^3 / 3 + x*t) dt
using Pkg;

Pkg.add("Plots")
using Plots;
n = 1.0;

function yfc(t)
    cos(1/3 * t^3 + n*t)
end;

# We're going to use the Runge-Kutta approximation, with the ODE:
# dx/dt = cos(t^3/3 + x*t)
# x(0) = 0
# t in 0 to tf (which should be high)
# This is a pretty poor approximation as the function doesn't converge, at least not well
const t0 = 0;
const tf = 1000;
const N = 1000000000;
const dt = 1000 / N;
t = t0:dt:tf;

x = 0;
for i in 1:N
    global x
    k1 = dt * yfc(t[i])
    k2 = dt * yfc(t[i] + 1/2*dt)
    k4 = dt * yfc(t[i+1])
    x = x + 1.0/6.0 * k1 + 2.0/3.0 * k2 + 1.0/6.0 * k4 
end
