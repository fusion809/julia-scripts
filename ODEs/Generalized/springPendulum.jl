using PyPlot;
using LaTeXStrings;
include("RKF45.jl");

# RHS of the problem expressed as a 1st order system
function springPen(params::NamedTuple, t::Float64, vars::SVector{4,Float64})::SVector{4,Float64}
    g = params.g;
    k = params.k;
    m = params.m;
    r = vars[1];
    dr = vars[2];
    th = vars[3];
    dth = vars[4];
    array = [dr, r*dth^2 - g*sin(th) - k*r/m, dth, -2*dr/r*dth - g/r*cos(th)];
    return array;
end

# Problem parameters
g = 9.81;
k = 10.0;
m = 1.0;
params = (g = g, k = k, m = m);

# Initial conditions
r0 = 1.0;
dr0 = 0.0;
theta0 = 0.0;
thetaDot0 = 0.0;
conds = @SVector [r0, dr0, theta0, thetaDot0];

# Error tolerance and step size
epsilon = 1e-8;
dtInitial = 0.1;

# Determine a suitable range for t
t0 = 0.0;
tf = 10.0;

# Solve problem
@time begin
t, vars = RKF45(springPen, params, t0, tf, conds, epsilon, dtInitial);
end
r = vars[:,1];
dr = vars[:,2];
th = vars[:,3];
dth = vars[:,4];
x = r.*cos.(th);
y = r.*sin.(th);

# Plot solution
PyPlot.figure(1);
PyPlot.clf();
PyPlot.plot(t, r, label=L"r");
PyPlot.legend();
PyPlot.figure(2);
PyPlot.clf();
PyPlot.plot(t, th, label=L"\theta");
PyPlot.legend();
PyPlot.figure(3);
PyPlot.clf();
PyPlot.plot(t, dr, label=L"\dot{r}");
# PyPlot.plot(t, dth, label=L"\dot{\theta}")
PyPlot.legend();
PyPlot.figure(4);
PyPlot.clf();
PyPlot.plot(t, dth, label=L"\dot{\theta}")
PyPlot.legend();
PyPlot.figure(5);
PyPlot.clf();
PyPlot.plot(x, y);
PyPlot.xlabel(L"x");
PyPlot.ylabel(L"y");
PyPlot.figure(6);
PyPlot.clf();
PyPlot.plot(th, dth)
PyPlot.xlabel(L"\theta")
PyPlot.ylabel(L"\dot{\theta}")
PyPlot.figure(7);
PyPlot.clf();
PyPlot.plot(r, dr)
PyPlot.xlabel(L"r");
PyPlot.ylabel(L"\dot{r}")