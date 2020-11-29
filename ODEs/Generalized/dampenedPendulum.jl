using PyPlot;
using LaTeXStrings;
include("RKF45.jl");

# RHS of the problem expressed as a 1st order system
function SP(params::NamedTuple, t::Float64, vars::SVector{2,Float64})::SVector{2,Float64}
    g = params.g;
    l = params.l;
    mu = params.mu;
    x = vars[1];
    dx = vars[2];
    array = [dx, -mu*dx-g/l*cos(x)];
    return array;
end

# Problem parameters
g = 9.81;
l = 1.0;
mu = 0.5;
params = (g = g, l = l, mu = mu);

# Initial conditions
theta0 = 0.0;
thetaDot0 = 10;
conds = @SVector [theta0, thetaDot0];

# Error tolerance and step size
epsilon = 3e-12;
dtInitial = 0.1;

# Determine a suitable range for t
t0 = 0.0;
tf = 10.0;

# Solve problem
@time begin
t, vars = RKF45(SP, params, t0, tf, conds, epsilon, dtInitial);
end
theta = vars[:,1];
thetaCor = rem.(theta, 2*pi);
thetaDot = vars[:,2];

# Plot solution
PyPlot.figure(1);
PyPlot.clf();
PyPlot.plot(t, thetaCor, label=L"\theta corrected");
# PyPlot.plot(t, theta, label=L"\theta");
PyPlot.plot(t, thetaDot, label=L"\dot{\theta}")
PyPlot.legend();
PyPlot.figure(2);
PyPlot.clf();
PyPlot.plot(thetaCor, thetaDot);
PyPlot.xlabel(L"\theta");
PyPlot.ylabel(L"\dot{\theta}");