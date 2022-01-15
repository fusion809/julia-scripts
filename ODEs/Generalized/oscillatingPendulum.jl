using PyPlot;
using LaTeXStrings;
include("RKF45.jl");

# RHS of the problem expressed as a 1st order system
function OscPen(params::NamedTuple, t::Float64, vars::SVector{4,Float64})::SVector{4,Float64}
    g = params.g;
    r0 = params.r0;
    k = params.k;
    m = params.m;
    r = vars[1];
    dr = vars[2];
    theta = vars[3];
    dtheta = vars[4];
    array = [dr, r*dtheta^2 - g*cos(theta) - k/m*(r-r0), dtheta, -2*(dr/r)*dtheta + g/r*sin(theta)];
    return array;
end

# Problem parameters
g = 9.81;
r0 = 1.0;
k = 1;
m = 1;
t0 = 0.0;
tf = 100.0;
params = (g = g, r0 = r0, k = k, m = m);

# Initial conditions
theta0 = pi/2;
thetaDot0 = 0.0;
dr0 = 1.0;
conds = @SVector [r0, dr0, theta0, thetaDot0];

# Error tolerance and step size
epsilon = 1e-10;
dtInitial = 0.1;

# Solve problem
t, vars = RKF45(OscPen, params, t0, tf, conds, epsilon, dtInitial);
r = vars[:,1];
dr = vars[:,2];
theta = vars[:,3];
thetaDot = vars[:,4];

# Plot solution
PyPlot.figure(1);
PyPlot.clf();
PyPlot.plot(t, vars[:,1], label=L"$r$");
PyPlot.plot(t, vars[:,2], label=L"$\frac{dr}{dt}$");
PyPlot.xlabel(L"$t$")
PyPlot.legend();
PyPlot.figure(2);
PyPlot.clf();
PyPlot.plot(t, vars[:,3], label=L"$\theta$");
PyPlot.plot(t, vars[:,4], label=L"$\frac{d\theta}{dt}$");
PyPlot.xlabel(L"$t$")
PyPlot.legend();
PyPlot.figure(3);
PyPlot.clf();
PyPlot.plot(vars[:,1], vars[:,2]);
PyPlot.xlabel(L"$r$")
PyPlot.ylabel(L"$\frac{dr}{dt}$", rotation=0)
PyPlot.figure(4);
PyPlot.clf();
PyPlot.plot(vars[:,3], vars[:,4]);
PyPlot.xlabel(L"$\theta$")
PyPlot.ylabel(L"$\frac{d\theta}{dt}$", rotation=0)