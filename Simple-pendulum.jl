# This script solves the problem of the double pendulum, but does not plot it
# It exists, in case I wish to try out other plotting methods than PyPlot
# Import package manager
using Pkg;

# Install and import ODE
Pkg.add("ODE")
using ODE;

function f(t, r)
	# Extract the coordinates from the r vector
	(x, dx) = r ;

	# The double pendulum equations
	dx_dt   = dx ;
	d2x_dt2 = - g/l * cos(x);
	# Return the derivatives as a vector
	[dx_dt; d2x_dt2]
end;

# Constants, I do like that I do not have to parse them manually to ode78
const g     = 9.8   # Acceleration due to gravity in m/s^2
const l     = 1     # Length of pendulum 1 in metres

# Define time vector and interval grid
const dt    = 1e-5       # Time increments
const tf    = 1e2        # Time at which the simulation ends
# A column vector containing discrete, dt-spaced, time values between 0 and tf
t           = 0:dt:tf

# Initial position in phase space
const r0    = [0.0; 0.0]

# Solve the equations
(t, pos)    = ode78(f, r0, t)
# N cannot be defined before ode78, as it redefines it
N           = length(t)

# Map pos values to x, dx, y and dy
x           = map(v -> v[1], pos)
dx          = map(v -> v[2], pos)
# Plots produces non-interative plots, unlike PyPlot
Pkg.add("PyPlot")
Pkg.add("LaTeXStrings")

using PyPlot;
using LaTeXStrings
PyPlot.figure(1)
PyPlot.plot(t,x)
PyPlot.xlabel(L"$t$")
PyPlot.ylabel(L"$\theta$")
PyPlot.figure(2)
PyPlot.plot(t,dx)
PyPlot.xlabel(L"$t$")
PyPlot.ylabel(L"$\dot{\theta}$")
PyPlot.figure(3)
PyPlot.plot(x,dx)
PyPlot.xlabel(L"$\theta$")
PyPlot.ylabel(L"$\dot{\theta}$")
PyPlot.figure(4)
PyPlot.plot3D(t,x,dx)
PyPlot.xlabel(L"$t$")
PyPlot.ylabel(L"$\theta$")
PyPlot.zlabel(L"$\dot{\theta}$")

# Write t, x, dx, y and dy to a file
A           = hcat(t, x, dx);
#Af          = A[N,:];
open("/data/GitHub/mine/scripts/julia-scripts/Simple-pendulum-dt-$dt-tf-$tf.txt", "w") do file
	write(file, "$A");
end
