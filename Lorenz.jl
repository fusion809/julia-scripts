# This is largely copied from the ODE.jl repository
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
const dt = 0.001
const tf = 100.0
t = 0:dt:tf

# Initial position in space
const r0 = [-2.0; 3.0; 0.5]

# Constants sigma, rho and beta
const sigma = 10.0
const rho   = 28.0
const bet   = 8.0/3.0;

(t, pos) = ode78(f, r0, t)
x = map(v -> v[1], pos)
y = map(v -> v[2], pos)
z = map(v -> v[3], pos);

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
