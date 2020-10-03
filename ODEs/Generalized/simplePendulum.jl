using FunctionIntegrator;
using PyPlot;
include("RKF45.jl");

# RHS of the problem expressed as a 1st order system
function f(params, t::Float64, vars::Vector)
    g = params.g;
    l = params.l;
    x = vars[1];
    dx = vars[2];
    array = [dx, -g/l*cos(x)];
    return array;
end

# Parameter object
struct paramObj
    g::Float64
    l::Float64
end

# Problem parameters
g = 9.81;
l = 1.0;
params = paramObj(g, l);

# Initial conditions
theta0 = 1.570796;
thetaDot0 = 0;

# Error tolerance and step size
epsilon = 3e-12;
dtInitial = 0.1;

# Determine a suitable range for t
if ( thetaDot0^2 + 2*g/l * (sin(theta0)-1) > 0)
    t0 = 0.0;
    tf = 10.0;
else
    t0 = 0.0;
    thetaMin = asin((thetaDot0^2+2*(g/l)*sin(theta0))/(2*g/l));
    thetaMax = - thetaMin - pi;
    tf = t0 + 4*abs(chebyshev_quadrature(theta -> 1/sqrt(thetaDot0^2+(2*g/l)*(sin(theta0)-sin(theta))), 1e5, 1, thetaMin, thetaMax));
end

# Solve problem
solution = RKF45(f, params, t0, tf, [theta0 thetaDot0], epsilon, dtInitial);

# Plot solution
PyPlot.figure(1);
PyPlot.plot(solution.t, solution.vars);
PyPlot.figure(2);
PyPlot.plot(solution.vars[:,1], solution.vars[:,2]);