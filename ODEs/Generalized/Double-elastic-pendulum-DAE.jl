using DifferentialEquations
using DAEProblemLibrary

function elastic_double_pendulum!(du, u, p, t)
    # Unpack variables
    r1, dr1, th1, dth1, r2, dr2, th2, dth2 = u
    m1, m2, g, k1, k2, l1, l2, b1, c1, b2, c2 = p
    Δ = th2 - th1
    v1  = sqrt(dr1^2+r1^2*dth1^2);
    v2  = sqrt(v1^2 + dr2^2 + r2^2*dth2^2 + 2*cos(Δ)*(dr1*dr2 + r1*r2*dth1*dth2)+2*sin(Δ)*(r1*dr2*dth1-dr1*r2*dth2));
    Qr1 = -(b1 + c1*v1)*dr1 - (b2 + c2*v2)*(dr1 + dr2*cos(Δ)-r2*dth2*sin(Δ));
    Qr2 = -(b2 + c2*v2)*(dr1*cos(Δ) + r1*dth1*sin(Δ)+dr2);
    Qt1 = -(b1 + c1*v1)*r1^2*dth1 - (b2 + c2*v2)*(r1^2*dth1 + r1*dr2*sin(Δ) + r1*r2*dth2*cos(Δ))
    Qt2 = -(b2 + c2*v2)*(r2^2*dth2 - dr1*r2*sin(Δ) + r1*r2*dth1*cos(Δ))

    # Position derivatives
    du[1] = dr1
    du[3] = dr2
    du[5] = dth1
    du[7] = dth2

    # Force equations (from your full system)
    # Use the original non-rearranged equations here, e.g.:

    # Equation for d²r₁
    du[2] = r1*dth1^2 - g*sin(th1) +
            (m2/(m1 + m2)) * (cos(Δ)*(r2*dth2^2 - du[4]) + sin(Δ)*(2*dr2*dth2 + r2*du[8])) -
            (Qr1-k1*(r1 - l1))/(m1 + m2)

    # Equation for d²r₂
    du[4] = -cos(Δ)*(du[2] - r1*dth1^2) - sin(Δ)*(r1*du[6] + 2*dr1*dth1) + r2*dth2^2 -
            g*sin(th2) + (Qr2 - k2*(r2 - l2))/m2

    # Equation for d²θ₁
    du[6] = -2*dr1*dth1/r1 - g*cos(th1)/r1 -
            (m2/((m1 + m2)*r1))*(cos(Δ)*(2*dr2*dth2 + r2*du[8]) + sin(Δ)*(du[4] - r2*dth2^2)) +
            Qt1/((m1 + m2)*r1^2)

    # Equation for d²θ₂
    du[8] = -2*dr2*dth2/r2 - (cos(Δ)/r2)*(2*dr1*dth1 + r1*du[6]) -
            (sin(Δ)/r2)*(r1*dth1^2 - du[2]) - g*cos(th2)/r2 + Qt2/(m2*r2^2)
end

# Initial conditions
u0 = [1, 0, 1, 0, 0, 0, 0, 0]
du0 = zeros(8)  # du[2], du[4], du[6], du[8] can be guessed/estimated

# Parameters: [m1, m2, g, k1, k2, l1, l2]
p = [1.0, 1.0, 9.81, 10.0, 10.0, 1.0, 1.0, 0.1, 0.02, 0.1, 0.02]

tspan = (0.0, 1.0)

prob = DAEProblem(elastic_double_pendulum!, du0, u0, tspan, p)
using Sundials;
sol = solve(prob, IDA())