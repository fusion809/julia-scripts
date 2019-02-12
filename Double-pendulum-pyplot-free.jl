# This script solves the problem of the double pendulum, but does not plot it
# It exists, in case I wish to try out other plotting methods than PyPlot
# Import package manager
using Pkg;

# Install and import ODE
Pkg.add("ODE")
using ODE;

function f(t, r)
	# Extract the coordinates from the r vector
	(x, dx, y, dy) = r

	# The double pendulum equations
	dx_dt   = dx
	d2x_dt2 = - ((g * (2 * m1 + m2) * sin(x) + m2 * (g * sin(x-2*y) + 2*(l2 * (dy^2) + l1 * (dx^2) * cos(x-y)) * sin(x-y)))/(2 *l1 * (m1 + m2 - m2 * cos(x-y)^2)))
	dy_dt   = dy
	d2y_dt2 = (((m1 + m2) * (l1 * (dx^2) + g * cos(x)) + l2*m2*(dy^2)*cos(x-y)) * sin(x-y)) / (l2 * (m1 + m2 - m2 * cos(x-y)^2));

	# Return the derivatives as a vector
	[dx_dt; d2x_dt2; dy_dt; d2y_dt2]
end;

# Constants, I do like that I do not have to parse them manually to ode78
const m1    = 1     # Mass of bob 1 in kg
const m2    = 1     # Mass of pendulum bob 2 in kg
const g     = 9.8   # Acceleration due to gravity in m/s^2
const l1    = 1     # Length of pendulum 1 in metres
const l2    = 1     # Length of pendulum 2 in metres

# Define time vector and interval grid
const dt    = 0.001    # Time increments
const tf    = 10.0     # Time at which the simulation ends
# A column vector containing discrete, dt-spaced, time values between 0 and tf
t           = 0:dt:tf

# Initial position in phase space
const r0    = [pi; 0.0; pi/2; 0.0]

# Solve the equations
(t, pos)    = ode78(f, r0, t)
# N cannot be defined before ode78, as it redefines it
N           = length(t)

# Map pos values to x, dx, y and dy
x           = map(v -> v[1], pos)
xf          = x[end];
dx          = map(v -> v[2], pos)
dxf         = dx[end];
y           = map(v -> v[3], pos)
yf          = y[end];
dy          = map(v -> v[4], pos);
dyf         = dy[end];
# Plots produces non-interative plots, unlike PyPlot
Pkg.add("Plots")
using Plots;
plot(x,y)

# Write t, x, dx, y and dy to a file
A           = hcat(t, x, dx, y, dy);
#Af          = A[N,:];
open("/data/GitHub/mine/math/julia-scripts/Double-pendulum-dt-$dt-tf-$tf.txt", "w") do file
	write(file, "$A");
end
