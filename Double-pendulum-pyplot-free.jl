# This script solves the problem of the double pendulum, but does not plot it
# It exists, in case I wish to try out other plotting methods than PyPlot
using ODE;

function f(t, r)
    # Extract the coordinates from the r vector
    (x, dx, y, dy) = r

    # The Lorenz equations
    dx_dt  = dx
    d2x_dt = - ((g * (2 * m1 + m2) * sin(x) + m2 * (g * sin(x-2*y) + 2*(l2 * (dy^2) + l1 * (dx^2) * cos(x-y)) * sin(x-y)))/(2 *l1 * (m1 + m2 - m2 * cos(x-y)^2)))
    dy_dt  = dy
    d2y_dt = (((m1 + m2) * (l1 * (dx^2) + g * cos(x)) + l2*m2*(dy^2)*cos(x-y)) * sin(x-y)) / (l2 * (m1 + m2 - m2 * cos(x-y)^2));

    # Return the derivatives as a vector
    [dx_dt; d2x_dt; dy_dt; d2y_dt]
end;

const m1   = 1
const m2   = 1
const g    = 9.8
const l1   = 1
const l2   = 1

# Define time vector and interval grid
const dt   = 0.001
const tf   = 1000.0
t          = 0:dt:tf

# Initial position in space
const r0   = [pi; 0.0; pi/2; 0.0]

(t, pos)   = ode78(f, r0, t)
x          = map(v -> v[1], pos)
dx         = map(v -> v[2], pos)
y          = map(v -> v[3], pos)
dy         = map(v -> v[4], pos);