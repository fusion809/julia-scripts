# This uses the equation:
# Ai(x) = int_0^inf cos(t^3 / 3 + x*t) dt
using Pkg;
Pkg.add("SpecialFunctions");
using SpecialFunctions;
x = 1.0;

function yfc(x,t)
    cos(1.0/3.0 * t^3 + x*t)/pi
end

function Simpson(dt,x,t,i,N)
    if i==1 || i == N+1
        return dt/3.0*yfc(x,t)
    elseif (i % 2) == 1
        return 2.0*dt/3.0*yfc(x,t)
    else
        return 4.0*dt/3.0*yfc(x,t)
    end
end
# We're going to use the Runge-Kutta approximation, with the ODE:
# dx/dt = cos(t^3/3 + x*t)
# x(0) = 0
# t in 0 to tf (which should be high)
# This is a pretty poor approximation as the function doesn't converge, at least not well
t0 = 0.0;
tf = 10000.0;
N = 100000000000;
dt = (tf-t0) / N;
t = t0:dt:tf;

function integrator(dt,x,t,N)
    y = 0;
    for i in 1:N+1
        y = y + Simpson(dt,x,t[i],i,N)
    end

    return y
end

y = integrator(dt,x,t,N);

error = abs(y-airyai(x));