using FunctionIntegrator;
using PyPlot;
include("RKF45.jl");

# RHS of the problem expressed as a 1st order system
function SP(params::NamedTuple, t::Float64, vars::Array{Float64,1})
    g = params.g;
    l = params.l;
    x = vars[1];
    dx = vars[2];
    array = [dx, -g/l*cos(x)];
    return array;
end

# Problem parameters
g = 9.81;
l = 1.0;
params = (g = g, l = l);

# Initial conditions
theta0 = 0.0;
thetaDot0 = sqrt(2*g/l)-1e-10;
conds = [theta0; thetaDot0];

# Error tolerance and step size
epsilon = 3e-12;
dtInitial = 0.1;

# Determine a suitable range for t
if ( thetaDot0^2 + 2*g/l * (sin(theta0)-1) > 0) # Condition for non-periodicity
    t0 = 0.0;
    tf = 10.0;
else
    t0 = 0.0;
    # Time when thetaDot = 0 that is closest to t0 (potentially in the past)
    thetaMin = asin((thetaDot0^2+2*(g/l)*sin(theta0))/(2*g/l));
    # Time when thetaDot = 0 next
    thetaMax = - thetaMin - pi;
    tf = t0 + 4*abs(chebyshev_quadrature(theta -> 1/sqrt(thetaDot0^2+(2*g/l)*(sin(theta0)-sin(theta))), 1e5, 1, thetaMin, thetaMax));
end

# Solve problem
solution = RKF45(SP, params, t0, tf, conds, epsilon, dtInitial);
vars = solution.vars;
t = solution.t;
theta = vars[:,1];
thetaDot = vars[:,2];

# Plot solution
PyPlot.figure(1);
PyPlot.clf();
PyPlot.plot(t, vars);
PyPlot.figure(2);
PyPlot.clf();
PyPlot.plot(vars[:,1], vars[:,2]);