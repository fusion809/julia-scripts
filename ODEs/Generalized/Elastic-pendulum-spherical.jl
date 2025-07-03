include("RKF45.jl")

# Euler-Lagrange equations
function ELE(params::NamedTuple, t::Float64, coords::SVector{6, Float64})::SVector{6,Float64}
    l          = params.l;
    g          = params.g;
    m          = params.m;
    k          = params.k;
    b1         = params.b1;
    c1         = params.c1;
    threshold  = params.th;
    kd         = k/m;
    xi         = coords[1];
    xidot      = coords[2];
    theta      = coords[3];
    thetadot   = coords[4];
    varphi     = coords[5];
    varphidot  = coords[6];
    r          = l + xi;
    v          = sqrt(xidot^2 + r^2*(thetadot^2*sin(varphi)^2+varphidot^2));
    Qxi        = -(b1+c1*v)*xidot;
    Qtheta     = -(b1+c1*v)*r^2*thetadot*(sin(varphi))^2;
    Qvarphi    = -(b1+c1*v)*r^2*varphidot;
    xiddot     = r*(thetadot^2*(sin(varphi))^2+varphidot^2) - kd*xi - g*cos(varphi) + Qxi/m;
    
    if (abs(sin(varphi)) < threshold)
        thetaddot = -2*thetadot * xidot/r + Qtheta/(m*r^2);
    else
        thetaddot = -2*thetadot*(xidot/r + varphidot*cot(varphi)) + Qtheta/(m*r^2*(sin(varphi))^2)
    end
    varphiddot = (g*sin(varphi)-2*xidot*varphidot)/r + thetadot^2*sin(2*varphi)/2 + Qvarphi/(m*r^2);
    return [xidot, xiddot, thetadot, thetaddot, varphidot, varphiddot]
end

# Hamilton's equations
function HE(params::NamedTuple, t::Float64, phaseSpace::SVector{6, Float64})::SVector{6, Float64}
    l            = params.l;
    g            = params.g;
    m            = params.m;
    k            = params.k;
    b1           = params.b1;
    c1           = params.c1;
    threshold    = params.th;
    kd           = k/m;
    xi           = phaseSpace[1];
    pxi          = phaseSpace[2];
    theta        = phaseSpace[3];
    ptheta       = phaseSpace[4];
    varphi       = phaseSpace[5];
    pvarphi      = phaseSpace[6];
    xidot        = pxi/m;
    r            = l+xi;
    if (abs(sin(varphi)) < threshold)
        thetadot = ptheta/(m*r^2)
    else
        thetadot = ptheta/(m*r^2*(sin(varphi))^2)
    end
    varphidot    = pvarphi/(m*r^2);
    if (abs(sin(varphi)) < threshold)
        pxidot   = ptheta^2/(m*r^3)
    else
        pxidot   = ptheta^2/(m*r^3*(sin(varphi))^2)
    end
    v            = sqrt(pxi^2/m^2+ptheta^2/(m^2*r^2*(sin(varphi)^2))+pvarphi^2/(m^2*r^2))
    Qxi          = -(b1+c1*v)*xidot;
    Qtheta       = -(b1+c1*v)*r^2*thetadot*(sin(varphi))^2;
    Qvarphi      = -(b1+c1*v)*r^2*varphidot;
    pxidot      += pvarphi^2/(m*r^3)  - k*xi - m*g*cos(varphi) + Qxi;
    pthetadot    = Qtheta;
    if (abs(sin(varphi)) < threshold)
        pvarphidot = ptheta^2/(m*r^2)
    else
        pvarphidot = ptheta^2*cot(varphi)*(csc(varphi))^2/(m*r^2)
    end
    pvarphidot += m*g*r*sin(varphi) + Qvarphi;
    return [xidot, pxidot, thetadot, pthetadot, varphidot, pvarphidot];
end

params         = (m=1, g=9.81, l=1, k=1, b1=0.0, c1=0.00, th=1e-5)
t0             = 0.0;
tf             = 300.0;
xi0            = 0.0;
xidot0         = 0.0;
r0             = params.l+xi0;
pxi0           = params.m*xidot0;
theta0         = 0.0;
thetadot0      = 0.0;
varphi0        = pi/2;
varphidot0     = 0.0;
ptheta0        = params.m*thetadot0*r0^2*sin(varphi0)^2;
pvarphi0       = params.m * varphidot0 * r0^2;
conds          = @SVector [xi0, xidot0, theta0, thetadot0, varphi0, varphidot0]
qpconds        = @SVector [xi0, pxi0, theta0, ptheta0, varphi0, pvarphi0];
epsilon        = 1e-10;
dtInit         = 1e-3;
t1, coords     = RKF45(ELE, params, t0, tf, conds, epsilon, dtInit, 
tolType="relative", dtMin=1e-3)
t2, phaseSpace = RKF45(HE, params, t0, tf, qpconds, epsilon, dtInit, 
tolType="relative", dtMin=1e-3)
using PyPlot
if @isdefined(p1)
    PyPlot.close(p1)
end
p1 = PyPlot.figure(1)
PyPlot.plot(t1, coords[:,1])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\xi", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\xi\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum ELE xi vs t plot.png")
if @isdefined(p2)
    PyPlot.close(p2)
end
p2 = PyPlot.figure(2)
PyPlot.plot(t1, coords[:,2])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\dot{\xi}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\dot{\\xi}\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum ELE xidot vs t plot.png")
if @isdefined(p3)
    PyPlot.close(p3)
end
p3 = PyPlot.figure(3)
PyPlot.plot(coords[:,1], coords[:,2])
PyPlot.xlabel(L"\xi", fontsize=14)
PyPlot.ylabel(L"\dot{\xi}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\dot{\\xi}\$ vs \$\\xi\$.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum ELE xidot vs xi plot.png")
if @isdefined(p4)
    PyPlot.close(p4)
end
p4 = PyPlot.figure(4)
PyPlot.plot(coords[:,3], coords[:,4])
PyPlot.xlabel(L"\theta", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\dot{\\theta}\$ vs \$\\theta\$.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum ELE thetadot vs theta plot.png")
if @isdefined(p5)
    PyPlot.close(p5)
end
p5 = PyPlot.figure(5)
PyPlot.plot(coords[:,5], coords[:,6])
PyPlot.xlabel(L"\varphi", fontsize=14)
PyPlot.ylabel(L"\dot{\varphi}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\dot{\\varphi}\$ vs \$\\varphi\$.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum ELE varphidot vs varphi plot.png")
if @isdefined(p6)
    PyPlot.close(p6)
end
p6 = PyPlot.figure(6)
PyPlot.plot(t1, coords[:,3])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\theta", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\theta\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum ELE theta vs t plot.png")
if @isdefined(p7)
    PyPlot.close(p7)
end
p7 = PyPlot.figure(7)
PyPlot.plot(t1, coords[:,4])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\dot{\\theta}\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum ELE thetadot vs t plot.png")
if @isdefined(p8)
    PyPlot.close(p8)
end
p8 = PyPlot.figure(8)
PyPlot.plot(t1, coords[:,5])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\varphi", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\varphi\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum ELE varphi vs t plot.png")
if @isdefined(p9)
    PyPlot.close(p9)
end
p9 = PyPlot.figure(9)
PyPlot.plot(t1, coords[:,6])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\dot{\varphi}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\dot{\\varphi}\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum ELE varphidot vs t plot.png")
if @isdefined(p10)
    PyPlot.close(p10)
end
p10 = PyPlot.figure(10)
PyPlot.plot(t2, phaseSpace[:,1])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\xi", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\xi\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum Hamiltonian xi vs t plot.png")
if @isdefined(p11)
    PyPlot.close(p11)
end
p11 = PyPlot.figure(11)
PyPlot.plot(t2, phaseSpace[:,2])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"p_{\xi}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$p_{\\xi}\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum Hamiltonian pxi vs t plot.png")
if @isdefined(p12)
    PyPlot.close(p12)
end
p12 = PyPlot.figure(12)
PyPlot.plot(t2, phaseSpace[:,3])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\theta", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\theta\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum Hamiltonian theta vs t plot.png")
if @isdefined(p13)
    PyPlot.close(p13)
end
p13 = PyPlot.figure(13)
PyPlot.plot(t2, phaseSpace[:,4])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"p_{\theta}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$p_{\\theta}\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum Hamiltonian ptheta vs t plot.png")
if @isdefined(p14)
    PyPlot.close(p14)
end
p14 = PyPlot.figure(14)
PyPlot.plot(t2, phaseSpace[:,5])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\varphi", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$\\varphi\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum Hamiltonian varphi vs t plot.png")
if @isdefined(p15)
    PyPlot.close(p15)
end
p15 = PyPlot.figure(15)
PyPlot.plot(t2, phaseSpace[:,6])
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"p_{\varphi}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$p_{\\varphi}\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum Hamiltonian pvarphi vs t plot.png")
if @isdefined(p16)
    PyPlot.close(p16)
end
p16 = PyPlot.figure(16)
PyPlot.plot(phaseSpace[:,1], phaseSpace[:,2])
PyPlot.xlabel(L"\xi", fontsize=14)
PyPlot.ylabel(L"p_{\xi}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$p_{\\varphi}\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum Hamiltonian pxi vs xi plot.png")
if @isdefined(p17)
    PyPlot.close(p17)
end
p17 = PyPlot.figure(17)
PyPlot.plot(phaseSpace[:,3], phaseSpace[:,4])
PyPlot.xlabel(L"\theta", fontsize=14)
PyPlot.ylabel(L"p_{\theta}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$p_{\\theta}\$ vs \$\\theta\$.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum Hamiltonian ptheta vs theta plot.png")
if @isdefined(p18)
    PyPlot.close(p18)
end
p18 = PyPlot.figure(18)
PyPlot.plot(phaseSpace[:,5], phaseSpace[:,6])
PyPlot.xlabel(L"\varphi", fontsize=14)
PyPlot.ylabel(L"p_{\varphi}", rotation=0, fontsize=14)
PyPlot.title("Elastic pendulum: \$p_{\\varphi}\$ vs \$\\varphi\$.", fontsize=16)
PyPlot.savefig("graphics/Elastic pendulum Hamiltonian pvarphi vs varphi plot.png")

using CairoMakie

function spherical_to_cartesian(coords, l)
    xi     = coords[:, 1]
    theta  = coords[:, 3]
    varphi = coords[:, 5]
    r      = l .+ xi

    x = r .* sin.(varphi) .* cos.(theta)
    y = r .* sin.(varphi) .* sin.(theta)
    z = r .* cos.(varphi)

    return x, y, z
end

using Dierckx
function animate_elastic_pendulum(coords, tvec, l, N; savepath = "graphics/elastic_pendulum_3d.mp4")
    x, y, z = spherical_to_cartesian(coords, l)
    tuni    = LinRange(tvec[1], tvec[end], N)
    dt      = tuni[2]-tuni[1]
    splx    = Spline1D(tvec, x)
    sply    = Spline1D(tvec, y)
    splz    = Spline1D(tvec, z)
    xuni    = evaluate(splx, tuni)
    yuni    = evaluate(sply, tuni)
    zuni    = evaluate(splz, tuni)

    fig = CairoMakie.Figure(resolution = (800, 600))
    padding = 0.01
    ax = Axis3(fig[1, 1], xlabel = "x", ylabel = "y", zlabel = "z")
    limits!(ax, minimum(xuni)-padding, maximum(xuni)+padding, minimum(yuni)-padding, maximum(yuni)+padding, minimum(zuni)-padding, maximum(zuni)+padding)

    bob_pos = Point3f(xuni[1], yuni[1], zuni[1])
    bob_plot = scatter!(ax, [bob_pos], color = :red, markersize = 15)

    rod_plot = lines!(ax, Point3f[(0, 0, 0), bob_pos], color = :black)

    record(fig, savepath, 1:N; framerate = round(Int, 1/dt)) do i
        pos = Point3f(xuni[i], yuni[i], zuni[i])
        bob_plot[1] = [pos]
        rod_plot[1] = Point3f[(0, 0, 0), pos]
    end

    println("Saved animation to: $savepath")
end

function animate_elastic_pendulum_2d(coords, t1, l, N)
    x1, y, z1 = spherical_to_cartesian(coords, l)
    tf = t1[end];
    t0 = t1[1];
    #t1 = t;
    #N=1000;
    dt = (tf-t0)/N;
    t_uni = t0:dt:tf;
    splx1 = Spline1D(t1, x1)
    splz1 = Spline1D(t1, z1)
    x1_uni = evaluate(splx1, t_uni)
    z1_uni = evaluate(splz1, t_uni)

    # Create a new figure and axis
    fig = CairoMakie.Figure(resolution = (900, 900))
    #resize_to_layout!(fig)
    ax = Axis(fig[1, 1],
        title = "Elastic Pendulum",
        xlabel = "x", ylabel = "z",
        aspect = 1.0
    )

    # Initialise objects as vectors of 2D points
    line_data = [Point2f(0.0, 0.0), Point2f(x1_uni[1], z1_uni[1])]
    bob_data = [Point2f(x1_uni[1], z1_uni[1])]

    pendulum_line = lines!(ax, line_data, linewidth=2)
    bob = scatter!(ax, bob_data, color=:red, markersize=10)

    # Output path
    output_path = "graphics/elastic_pendulum_3d_in_2d.mp4"
    padding = 0.5
    xrange = extrema(x1)
    zrange = extrema(z1)

    xlims!(ax, xrange[1] - padding, xrange[2] + padding)
    ylims!(ax, zrange[1] - padding, zrange[2] + padding)

    # Animation loop
    record(fig, output_path, 1:N; framerate=round(Int, 1/dt)) do i
        # Update line and bob by replacing the full vectors of points
        pendulum_line[1][] = [Point2f(0.0, 0.0), Point2f(x1_uni[i], z1_uni[i])]
        bob[1][] = [Point2f(x1_uni[i], z1_uni[i])]
    end
end
N = 10000;
anim3d = animate_elastic_pendulum(coords, t1, params.l, N)
anim2d = animate_elastic_pendulum_2d(coords, t1, params.l, N)