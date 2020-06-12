# This uses the equation:
# Ai(x) = int_0^inf cos(t^3 / 3 + x*t) dt
using Pkg;
Pkg.add("QuadGK")
using QuadGK;
n = 1.0;

function yfc(t)
    cos(1.0/3.0 * t^3 + n*t)
end

function Simpson(dt,t,i,N)
    if i==0 || i == N
        return dt/3.0*yfc(t)
    elseif (i % 2) == 0
        return 2.0*dt/3.0*yfc(t)
    else
        return 4.0*dt/3.0*yfc(t)
    end
end
# We're going to use the Runge-Kutta approximation, with the ODE:
# dx/dt = cos(t^3/3 + x*t)
# x(0) = 0
# t in 0 to tf (which should be high)
# This is a pretty poor approximation as the function doesn't converge, at least not well
t0 = 0.0;
tf = 10000.0;
N = 1000000000;
dt = (tf-t0) / N;
t = t0:dt:tf;

x = 0;
for i in 1:N+1
    global x
    x = x + Simpson(dt,t[i],i,N)
end

xd, err = quadgk(t -> cos((1.0/3.0)*t^3+n*t), t0, tf, rtol=1e-10);