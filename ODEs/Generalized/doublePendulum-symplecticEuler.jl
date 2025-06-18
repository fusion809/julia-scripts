using LinearAlgebra, PyPlot

# Physical parameters
const m1 = 1.0
const m2 = 1.0
const r1 = 1.0
const r2 = 1.0
const g = 9.81

# Time parameters
tmax = 20.0
h = 0.00001
N = round(Int, tmax / h)

# Arrays
theta1 = zeros(N+1)
theta2 = zeros(N+1)
p1 = zeros(N+1)
p2 = zeros(N+1)
t = h .* (0:N)

# Initial conditions
theta1[1] = 0  # measured from x-axis
theta2[1] = 0
p1[1] = 0.0
p2[1] = 0.0

# Mass matrix function
function M(theta1, theta2)
    Delta= theta2 - theta1
    M11 = (m1 + m2) * r1^2
    M12 = m2 * r1 * r2 * cos(Delta)
    M21 = M12
    M22 = m2 * r2^2
    return [M11 M12; M21 M22]
end

# Gradient of the potential energy
function dV_dtheta(theta1, theta2)
    dV1 = (m1 + m2) * g * r1 * cos(theta1)
    dV2 = m2 * g * r2 * cos(theta2)
    return [dV1; dV2]
end
# Symplectic Euler loop
for n in 1:N
    # Current state
    th1 = theta1[n]
    th2 = theta2[n]
    p = [p1[n]; p2[n]]

    # Compute gradient of potential
    dV = dV_dtheta(th1, th2)

    # Momentum update (explicit)
    p_new = p - h .* dV

    # Invert mass matrix to get angular velocities
    Minv = inv(M(th1, th2))
    dtheta = Minv * p_new

    # Position update (implicit in p)
    theta1[n+1] = th1 + h * dtheta[1]
    theta2[n+1] = th2 + h * dtheta[2]

    # Save momenta
    p1[n+1] = p_new[1]
    p2[n+1] = p_new[2]
end

# Plot angles
PyPlot.figure(1)
PyPlot.plot(t, theta1)
PyPlot.xlabel("Time (s)")
PyPlot.ylabel(L"\theta_1", rotation=0)
PyPlot.title("Double Pendulum: Symplectic Euler Method")

PyPlot.figure(2)
PyPlot.plot(t, theta2)
PyPlot.xlabel("Time (s)")
PyPlot.ylabel(L"\theta_2", rotation=0)
PyPlot.title("Double Pendulum: Symplectic Euler Method")

PyPlot.figure(3)
PyPlot.plot(theta1, theta2)
PyPlot.xlabel(L"\theta_1")
PyPlot.ylabel(L"\theta_2", rotation=0)
PyPlot.title("Double Pendulum: Symplectic Euler Method")