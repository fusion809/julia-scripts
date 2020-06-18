# This is largely copied from the ODE.jl repository
# Import package manager
using Pkg;

# Install and import ODE
Pkg.add("ODE")
using ODE;

function f(t, r)
	# Extract the coordinates from the r vector
	(x, y, z) = r

	# The Lorenz equations
	dx_dt = sigma*(y - x)
	dy_dt = x*(rho - z) - y
	dz_dt = x*y - bet*z

	# Return the derivatives as a vector
	[dx_dt; dy_dt; dz_dt]
end;

# Define time vector and interval grid
dt = 0.001
tf = 1e2
t = 0:dt:tf

# Initial position in space
r0 = [-2.0; 3.0; 0.5]

# Constants sigma, rho and beta
sigma = 10.0
rho   = 28.0
bet   = 8.0/3.0;

(t, pos) = ode78(f, r0, t)
x = map(v -> v[1], pos)
y = map(v -> v[2], pos)
z = map(v -> v[3], pos);

# Get PyPlot and load it
Pkg.add("PyPlot")
using PyPlot;
#Pkg.add("Plots")
#Pkg.add("Iterators")
#using Plots

PyPlot.figure(1)
PyPlot.plot3D(x, y, z);
PyPlot.xlabel("x")
PyPlot.ylabel("y")
PyPlot.zlabel("z")

PyPlot.figure(2)
PyPlot.plot3D(t, x, y);
PyPlot.xlabel("t")
PyPlot.ylabel("x")
PyPlot.zlabel("y")

PyPlot.figure(3)
PyPlot.plot3D(t, y, z);
PyPlot.xlabel("t")
PyPlot.ylabel("y")
PyPlot.zlabel("z")

PyPlot.figure(4)
PyPlot.plot3D(t, x, z);
PyPlot.xlabel("t")
PyPlot.ylabel("x")
PyPlot.zlabel("z")

fig, ax = subplots(1, 3, sharex=true, sharey=true, figsize=(16,8))

ax[1][:plot](x, y)
ax[1][:set_title]("X-Y cut")

ax[2][:plot](x, z)
ax[2][:set_title]("X-Z cut")

ax[3][:plot](y, z)
ax[3][:set_title]("Y-Z cut");

#anim = @animate for i=1:length(t)
#    plot3D(x[i],y[i],z[i])
#end
#gif(anim, "/tmp/anim_fps15.gif", fps = 15)
