# This is largely copied from the ODE.jl repository
using ODE;

function f(t, r)
    # Extract the coordinates from the r vector
    (x, dx) = r

    # The Lorenz equations
    dX = dx
    d2X = -g/l*cos(x) + f0*cos(omega*t)

    # Return the derivatives as a vector
    [dX; d2X]
end;

# Define time vector and interval grid
const dt = 0.001
const tf = 100.0
t = 0:dt:tf

# Initial position in space
const r0 = [0.1; 0.0]

# Constants alpha, beta, gamma, delta and omega
const g     = 9.8
const l     = 1
const f0    = 1
const omega = 1

(t, pos) = ode78(f, r0, t)
x = map(v -> v[1], pos)
dx = map(v -> v[2], pos)

using PyPlot
figure(1)
plot(x, dx);
title("Phase plot");
xlabel("x")
ylabel("dx")

figure(2)
plot(t, x);
xlabel("t")
ylabel("x")

figure(3)
plot(t,dx);
xlabel("t")
ylabel("dx")
