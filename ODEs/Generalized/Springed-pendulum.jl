include("RKF45.jl")
params = (m1 = 1.0, m2 = 1.0, l=1.0, g=9.81, r=1, k=1, 
b1=0.0, c1=0.0, b2=0.0, c2=0.0);
z0      = 0.0; 
dz0     = 0.0;
pz0     = 0.0;
theta0  = 0.0;
dtheta0 = 0.0;
ptheta0 = 0.0;
t0 = 0.0;
tf = 100.0;
phSpconds = @SVector [z0, pz0, theta0, ptheta0];
conds = @SVector [z0, dz0, theta0, dtheta0];
dtInit = 1e-3;
epsilon = 1e-9;

function Hsystem(params::NamedTuple, t::Float64, phaseSpace::SVector{4, Float64})::SVector{4,Float64}
    z           = phaseSpace[1];
    pz          = phaseSpace[2];
    theta       = phaseSpace[3];
    ptheta      = phaseSpace[4];
    m1          = params.m1;
    m2          = params.m2;
    g           = params.g;
    l           = params.l;
    r           = params.r;
    k           = params.k;
    b1          = params.b1;
    c1          = params.c1;
    b2          = params.b2;
    c2          = params.c2;
    den         = m1+m2*(sin(theta))^2;
    dzL         = (pz-ptheta*cos(theta)/r)/den;
    dthetaLden  = m2*r^2*den;
    dthetaL     = 1/dthetaLden * (ptheta*(m1+m2)-m2*r*pz*cos(theta));
    dz          = 1/den * ((m1+m2)*dzL + m2*r*dthetaL*cos(theta));
    dthetadth   = 1/(r^2*den^2)*(pz*r*sin(theta)*(m1+m2*(1+(cos(theta))^2))-ptheta*(m1+m2)*sin(2*theta));
    dzdth       = 1/den^2 * (ptheta*sin(theta)/r * (m1+m2*(1+cos(theta)^2)) - m2*sin(2*theta)*pz);
    dzdpz       = 1/den;
    dzdpth      = -cos(theta)/(r*den);
    dthdpz      = dzdpth;
    dthdpth     = (m1+m2)/(m2*r^2*den)
    dtheta      = (m1+m2)*(dz*dzdpth) + m2*(r^2*dthetaL*dthdpth + r*dzdpth*dthetaL *cos(theta)+r*dz * dthdpth*cos(theta));

    v2          = sqrt(r^2*dtheta^2+dz^2+2*r*dz*dtheta*cos(theta));
    Qz          = -(b1+c1*abs(dz))*dz - (b2+c2*v2)*(dz+r*dtheta*cos(theta));
    Qtheta      = -(b2+c2*v2)*(r^2*dtheta + r*dz*cos(theta));
    dptheta     = -m2*(r^2*dtheta * dthetadth + r*dz * dthetadth *cos(theta) + r*dzdth*dtheta*cos(theta)-r*dz*dtheta*sin(theta)+g*r*cos(theta)) + Qtheta;
    dpz         = -(m1+m2)*g-k*z+Qz;

    return [dz, dpz, dtheta, dptheta];
end

function Lsystem(params::NamedTuple, t::Float64, coords::SVector{4,Float64})::SVector{4,Float64}
    z           = coords[1];
    dz          = coords[2];
    theta       = coords[3];
    dtheta      = coords[4];
    m1          = params.m1;
    m2          = params.m2;
    M           = m1+m2;
    g           = params.g;
    l           = params.l;
    r           = params.r;
    k           = params.k;
    b1          = params.b1;
    c1          = params.c1;
    b2          = params.b2;
    c2          = params.c2;
    alpha       = 1/(m1+m2*(sin(theta))^2);
    v2          = sqrt(r^2*dtheta^2+dz^2+2*r*dz*dtheta*cos(theta));
    Qz          = -(b1+c1*abs(dz))*dz - (b2+c2*v2)*(dz+r*dtheta*cos(theta));
    Qtheta      = -(b2+c2*v2)*(r^2*dtheta + r*dz*cos(theta));
    d2theta     = -alpha*m2*sin(theta)*cos(theta)*dtheta^2 - alpha*cos(theta)*(Qz-k*z)/r +alpha*Qtheta/(m2*r^2);
    d2z         = -g - m2*r/M *(d2theta*cos(theta)-dtheta^2*sin(theta))+(Qz-k*z)/M;
    return [dz, d2z, dtheta, d2theta];
end

# t, phaseSpace = RKF45(Hsystem, params, t0, tf, phSpconds, epsilon, dtInit);
# z         = phaseSpace[:,1];
# pz        = phaseSpace[:,2];
# theta     = phaseSpace[:,3];
# ptheta    = phaseSpace[:,4];
t, coords   = RKF45(Lsystem, params, t0, tf, conds, epsilon, dtInit)
z           = coords[:,1];
dz          = coords[:,2];
theta       = coords[:,3];
dtheta      = coords[:,4];
using PyPlot
if @isdefined(p1)
    PyPlot.close(p1)
end
p1 = PyPlot.figure(1)
PyPlot.plot(t, theta)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\theta", rotation=0, fontsize=14)
PyPlot.title("Single pendulum attached to a spring: \$\\theta\$ vs t", fontsize=16)
if @isdefined(p2)
    PyPlot.close(p2)
end
p2 = PyPlot.figure(2)
PyPlot.plot(t, z)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"z", rotation=0, fontsize=14)
PyPlot.title("Single pendulum attached to a spring: \$z\$ vs t", fontsize=16)

if @isdefined(p3)
    PyPlot.close(p3)
end
p3 = PyPlot.figure(3)
PyPlot.plot(z, dz)
PyPlot.xlabel(L"z", fontsize=14)
PyPlot.ylabel(L"\dot{z}", rotation=0, fontsize=14)
PyPlot.title("Single pendulum attached to a spring: \$\\dot{z}\$ vs z", fontsize=16)

if @isdefined(p4)
    PyPlot.close(p4)
end
p4 = PyPlot.figure(4)
PyPlot.plot(theta, dtheta)
PyPlot.xlabel(L"\theta", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}", rotation=0, fontsize=14)
PyPlot.title("Single pendulum attached to a spring: \$\\dot{\\theta}\$ vs \$\\theta\$", fontsize=16)

if @isdefined(p5)
    PyPlot.close(p5)
end
p5 = PyPlot.figure(5)
PyPlot.plot(z, dz)
PyPlot.xlabel(L"z", fontsize=14)
PyPlot.ylabel(L"\theta", rotation=0, fontsize=14)
PyPlot.title("Single pendulum attached to a spring: \$\\theta\$ vs z", fontsize=16)

x1 = params.r*cos.(theta);
y1 = params.l .+ z .+ params.r*sin.(theta);
using Dierckx, CairoMakie

function animate(t1, x1, y1, N)
    tf = t1[end];
    t0 = t1[1];
    dt = (tf-t0)/N;
    t_uni = t0:dt:tf;
    splx1 = Spline1D(t1, x1)
    sply1 = Spline1D(t1, y1)
    x1_uni = evaluate(splx1, t_uni)
    y1_uni = evaluate(sply1, t_uni)

    # Create a new figure and axis
    fig = CairoMakie.Figure(resolution = (800, 1200))
    #resize_to_layout!(fig)
    ax = Axis(fig[1, 1],
        title = "Springed pendulum",
        xlabel = "x", ylabel = "y",
        aspect = 1.0
    )

    # Initialise objects as vectors of 2D points
    line_data = [Point2f(0.0, 0.0), Point2f(x1_uni[1], y1_uni[1])]
    bob_data = [Point2f(x1_uni[1], y1_uni[1])]

    pendulum_line = lines!(ax, line_data, linewidth=2)
    bob = scatter!(ax, bob_data, color=:red, markersize=10)

    # Output path
    output_path = "graphics/springed_pendulumLag.mp4"
    padding = 0.5
    xrange = extrema(x1)
    yrange = extrema(y1)

    xlims!(ax, xrange[1] - padding, xrange[2] + padding)
    ylims!(ax, yrange[1] - padding, yrange[2] + padding)

    # Animation loop
    record(fig, output_path, 1:N; framerate=round(Int, 1/dt)) do i
        # Update line and bob by replacing the full vectors of points
        pendulum_line[1][] = [Point2f(0.0, 0.0), Point2f(x1_uni[i], y1_uni[i])]
        bob[1][] = [Point2f(x1_uni[i], y1_uni[i])]
    end
end

#animate(t, x1, y1, round(Int, tf*100))