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
    dx_dt = y - a*x.^3 + b*x.^2 - z + I
	dy_dt = c - d*x.^3 - y
	dz_dt = r*(s*(x-xr)-z)

	# Return the derivatives as a vector
	[dx_dt; dy_dt; dz_dt]
end;

# Define time vector and interval grid
const dt = 0.0001
const tf = 100.0
t = 0:dt:tf

# Initial position in space
const r0 = [0.7; 0.0; 0.0]

# Constants sigma, rho and beta
const a  = 1.0;
const b  = 3.0;
const c  = 1.0;
const d  = 5.0;
const r  = 1e-3;
const s  = 4.0;
const xr = -1.6;
const I  = -10.0;

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

#PyPlot.figure(2)
#PyPlot.plot3D(t, x, y);
#PyPlot.xlabel("t")
#PyPlot.ylabel("x")
#PyPlot.zlabel("y")

#PyPlot.figure(3)
#PyPlot.plot3D(t, y, z);
#PyPlot.xlabel("t")
#PyPlot.ylabel("y")
#PyPlot.zlabel("z")

#PyPlot.figure(4)
#PyPlot.plot3D(t, x, z);
#PyPlot.xlabel("t")
#PyPlot.ylabel("x")
#PyPlot.zlabel("z")
