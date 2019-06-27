# This is largely copied from the ODE.jl repository
# Import package manager
using Pkg;

# Install and import ODE
Pkg.add("ODE")
using ODE;

function f(t, r)
	# Extract the coordinates from the r vector
	(x, y, z) = r

	# The Rossler attractor
	dx_dt = -y-z
	dy_dt = x+a*y
	dz_dt = b+z*(x-c)

	# Return the derivatives as a vector
	[dx_dt; dy_dt; dz_dt]
end;

# Define time vector and interval grid
const dt = 1e-4
const tf = 1e3
t = 0:dt:tf

# Initial position in space
const r0 = [-0.1; 0.5; -0.6]

# Constants
const a     = 0.2
const b     = 0.2
const c     = 5.7

(t, pos) = ode78(f, r0, t)
x = map(v -> v[1], pos)
y = map(v -> v[2], pos)
z = map(v -> v[3], pos);

# Get PyPlot and load it
Pkg.add("PyPlot")
using PyPlot

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