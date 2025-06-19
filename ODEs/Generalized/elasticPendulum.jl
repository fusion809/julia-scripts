include("RKF45.jl")
function EP(params::NamedTuple, t::Float64, phaseSpace::SVector{4,Float64})::SVector{4,Float64}
    k         = params.k;
    m         = params.m;
    g         = params.g;
    l         = params.l;
    x         = phaseSpace[1];
    px        = phaseSpace[2];
    theta     = phaseSpace[3];
    ptheta    = phaseSpace[4];
    xdot      = px./m;
    pxdot     = (ptheta.^2)./(m*(x+l).^3) .- k*x .- m*g*sin.(theta);
    thetadot  = ptheta./(m*(x+l).^2);
    pthetadot = -m*g*(x+l).*cos(theta);
    return [xdot, pxdot, thetadot, pthetadot];
end

function Hamiltonian(params::NamedTuple, phaseSpace)
    k         = params.k;
    m         = params.m;
    g         = params.g;
    l         = params.l;
    x         = phaseSpace[:,1];
    px        = phaseSpace[:,2];
    theta     = phaseSpace[:,3];
    ptheta    = phaseSpace[:,4];
    H         = (px.^2)/(2*m) .+ (ptheta.^2)./(2*m*(x.+l).^2) .+ k*(x.^2)/2 .+ m*g*(x.+l).*sin.(theta);
    return H;
end

include("RKF45.jl")
m             = 1;
g             = 9.81;
l             = 0.5;
k             = 1;
params        = (m=m, g=g, l=l, k=k);
t0            = 0.0;
tf            = 60.0;
epsilon       = 1e-10;
dtInitial     = 1e-3;
x0            = 0.0;
px0           = 0.0;
theta0        = 0.0;
ptheta0       = 0.0;
phaseSpace0   = @SVector [x0, px0, theta0, ptheta0];
t, phaseSpace = RKF45(EP, params, t0, tf, phaseSpace0, epsilon, dtInitial);
x             = phaseSpace[:,1];
px            = phaseSpace[:,2];
theta         = phaseSpace[:,3];
ptheta        = phaseSpace[:,4];

using PyPlot
PyPlot.figure(1)
PyPlot.plot(t, x)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"x", rotation=0)
PyPlot.figure(2)
PyPlot.plot(t, theta)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"\theta", rotation=0)
PyPlot.figure(3)
PyPlot.plot(t, px)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"p_x", rotation=0)
PyPlot.figure(4)
PyPlot.plot(t, ptheta)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"p_{\theta}", rotation=0)
PyPlot.figure(5)
PyPlot.plot(x, theta)
PyPlot.xlabel(L"x")
PyPlot.ylabel(L"\theta", rotation=0)
PyPlot.figure(6)
PyPlot.plot(px, ptheta)
PyPlot.xlabel(L"p_x")
PyPlot.ylabel(L"p_{\theta}", rotation=0)
PyPlot.figure(7)
PyPlot.plot(x, px)
PyPlot.xlabel(L"x")
PyPlot.ylabel(L"p_x", rotation=0)
PyPlot.figure(8)
PyPlot.plot(theta, ptheta)
PyPlot.xlabel(L"\theta")
PyPlot.ylabel(L"p_{\theta}", rotation=0)
H = Hamiltonian(params, phaseSpace);
PyPlot.figure(9)
PyPlot.plot(t, H)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"H", rotation=0)
x1 = (x.+l).*cos.(theta);
y1 = (x.+l).*sin.(theta);
using CairoMakie, Dierckx
function animate_EP(t1, x1, y1, N)
    tf = t1[end];
    t0 = t1[1];
    dt = (tf-t0)/N;
    t_uni = t0:dt:tf;
    splx1 = Spline1D(t1, x1)
    sply1 = Spline1D(t1, y1)
    x1_uni = evaluate(splx1, t_uni)
    y1_uni = evaluate(sply1, t_uni)

    # Create a new figure and axis
    fig = CairoMakie.Figure(resolution = (9000, 900))
    ax = Axis(fig[1, 1],
        title = "Elastic Pendulum",
        xlabel = "x", ylabel = "y",
        aspect = DataAspect()
    )

    # Initialise objects as vectors of 2D points
    line_data = [Point2f(0.0, 0.0), Point2f(x1_uni[1], y1_uni[1])]
    bob_data = [Point2f(x1_uni[1], y1_uni[1])]

    pendulum_line = lines!(ax, line_data, linewidth=2)
    bob = scatter!(ax, bob_data, color=:red, markersize=10)

    # Output path
    output_path = "graphics/elastic_pendulum.gif"
    xlims!(ax, minimum(x1) - 0.5, maximum(x1) + 0.5)
    ylims!(ax, minimum(y1) - 0.5, maximum(y1) + 0.5)

    # Animation loop
    record(fig, output_path, 1:N; framerate=round(Int, 1/dt)) do i
        # Update line and bob by replacing the full vectors of points
        pendulum_line[1][] = [Point2f(0.0, 0.0), Point2f(x1_uni[i], y1_uni[i])]
        bob[1][] = [Point2f(x1_uni[i], y1_uni[i])]
    end
end
N=1000;
animate_EP(t, x1, y1, N)