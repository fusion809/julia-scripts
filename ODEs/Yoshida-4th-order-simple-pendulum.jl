using PyPlot

# Physical parameters
const m = 1.0      # mass (kg)
const r = 1.0      # pendulum length (m)
const g = 9.81     # gravity (m/s²)

# Simulation parameters
const tmax = 60.0
const h = 0.0001
const N = Int(tmax / h)

# Allocate arrays
θ = zeros(N+1)
p = zeros(N+1)
t = h .* (0:N)

# Initial conditions
θ[1] = 0    # initial angle (from horizontal)
p[1] = 0.0      # initial angular momentum

# Yoshida 4th-order coefficients
w1 = 1 / (2 - 2^(1/3))
w2 = -2^(1/3) / (2 - 2^(1/3))
a = [w1/2, (w2 + w1)/2, (w2 + w1)/2, w1/2]
b = [w1, w2, w1]

# Derivatives
function V′(θ)
    return m * g * r * cos(θ)
end

function force(θ)
    return -V′(θ)  # = -mgr cos(θ)
end

# Integrator loop
for n in 1:N
    θn = θ[n]
    pn = p[n]

    for i in 1:3
        pn += a[i] * h * force(θn)
        θn += b[i] * h * pn / (m * r^2)
    end
    pn += a[4] * h * force(θn)

    θ[n+1] = θn
    p[n+1] = pn
end

# Plot angle vs time
PyPlot.figure(1)
PyPlot.plot(t, θ)
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("Angle (rad)")
PyPlot.title("Simple Pendulum: Yoshida 4th-Order Symplectic Integrator")
PyPlot.figure(2)
PyPlot.plot(θ, p)
PyPlot.xlabel("θ")
PyPlot.ylabel("p")
PyPlot.figure(3)
PyPlot.plot(t, p)
PyPlot.xlabel("t")
PyPlot.ylabel("p")
# Plot energy vs time (optional)
H = p.^2 ./ (2 * m * r^2) .+ m * g * r .* sin.(θ)
PyPlot.figure(4)
PyPlot.plot(t, H)
PyPlot.title("Hamiltonian over time")