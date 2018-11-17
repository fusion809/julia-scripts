# This script solves the problem of the double pendulum
using ODE;

function f(t, r)
	# Extract the coordinates from the r vector
	(x, dx, y, dy) = r

	# The double pendulum equations
	dx_dt  = dx
	d2x_dt = - ((g * (2 * m1 + m2) * sin(x) + m2 * (g * sin(x-2*y) + 2*(l2 * (dy^2) + l1 * (dx^2) * cos(x-y)) * sin(x-y)))/(2 *l1 * (m1 + m2 - m2 * cos(x-y)^2)))
	dy_dt  = dy
	d2y_dt = (((m1 + m2) * (l1 * (dx^2) + g * cos(x)) + l2*m2*(dy^2)*cos(x-y)) * sin(x-y)) / (l2 * (m1 + m2 - m2 * cos(x-y)^2));

	# Return the derivatives as a vector
	[dx_dt; d2x_dt; dy_dt; d2y_dt]
end;

# Constants, I do like that I do not have to parse them manually to ode78
const m1   = 1
const m2   = 1
const g    = 9.8
const l1   = 1
const l2   = 1

# Define time vector and interval grid
const dt   = 0.001
const tf   = 1000.0
# A column vector containing discrete, dt-spaced, time values between 0 and tf
t          = 0:dt:tf

# Initial position in phase space
const r0   = [pi; 0.0; pi/2; 0.0]

# Solve the equations
(t, pos)   = ode78(f, r0, t)

# Map pos values to x, dx, y and dy
x          = map(v -> v[1], pos)
dx         = map(v -> v[2], pos)
y          = map(v -> v[3], pos)
dy         = map(v -> v[4], pos);

# Plot the data, against time and in phase space
using PyPlot
PyPlot.figure(1)
PyPlot.plot(x, y)
PyPlot.xlabel(L"$\theta_1$")
PyPlot.ylabel(L"$\theta_2$")

PyPlot.figure(2)
PyPlot.plot(x, dx)
PyPlot.xlabel(L"$\theta_1$")
PyPlot.ylabel(L"$\dot{\theta_1}$")

PyPlot.figure(3)
PyPlot.plot(y, dy)
PyPlot.xlabel(L"$\theta_2$")
PyPlot.ylabel(L"$\dot{\theta_2}$")

PyPlot.figure(4)
PyPlot.plot(dx,dy)
PyPlot.xlabel(L"$\dot{\theta_{1}}$")
PyPlot.ylabel(L"$\dot{\theta_{2}}$")

PyPlot.figure(5)
PyPlot.plot(t, x)
PyPlot.xlabel("t")
PyPlot.ylabel(L"$\theta_1$")

PyPlot.figure(6)
PyPlot.plot(t, y)
PyPlot.xlabel("t")
PyPlot.ylabel(L"$\theta_2$")

PyPlot.figure(7)
PyPlot.plot3D(t, x, y)
PyPlot.xlabel("t")
PyPlot.ylabel("x")
PyPlot.zlabel("y")
