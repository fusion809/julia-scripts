# This is largely copied from the ODE.jl repository
# Import package manager
using Pkg;

# Install and import ODE
Pkg.add("ODE")
using ODE;

function f(t, r)
    # Extract the coordinates from the r vector
    (x, dx) = r

    # The Lorenz equations
    dX = dx
    d2X = gamma*cos(omega*t) - beta*x - delta*dx - alpha*x^3

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
const alpha = 1
const beta  = -1
const delta = 0.2
const gamma = 0.3
const omega = 1

(t, pos) = ode78(f, r0, t)
x = map(v -> v[1], pos)
dx = map(v -> v[2], pos)

# Get PyPlot and load it
Pkg.add("PyPlot")
using PyPlot

figure(1)
plot(x, dx);
title("Phase plot");
xlabel("x")
ylabel(L"$\dot{x}$")

figure(2)
plot(t, x);
title(L"$x$ vs $t$")
xlabel("t")
ylabel("x")

figure(3)
plot(t,dx);
title(L"$\dot{x}$ vs $t$")
xlabel("t")
ylabel(L"$\dot{x}$")

figure(4)
plot3D(t,x,dx)
title(L"Phase plot vs $t$")
xlabel("t")
ylabel("x")
zlabel(L"$\dot{x}$")
