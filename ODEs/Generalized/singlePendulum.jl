include("RKF45.jl")
params = (mr = 1, mb = 1, g=9.81, l = 1.0, br = 0.1, cr = 0.04, bb = 0.1, cb = 0.04)
theta0 = 0.0;
dtheta0 = 0.0;
conds = @SVector [theta0, dtheta0];
t0 = 0.0;
tf = 60.0;
epsilon = 1e-11;
dtInit = 1e-3;

function sinPen(params::NamedTuple, t::Float64, coords::SVector{2,Float64})::SVector{2,Float64}
    mr = params.mr;
    mb = params.mb;
    g = params.g;
    l = params.l;
    br = params.br;
    bb = params.bb;
    cr = params.cr;
    cb = params.cb;
    dtheta = coords[2];
    theta = coords[1];
    mv = mr+2*mb;
    mt = mr+3*mb;
    d2theta = -3*mv*g/(2*mt*l)*cos(theta) - 3*dtheta/(8*mt)*(2*br + 8*bb + l*dtheta*(cr + 8*cb));
    return [dtheta, d2theta];
end

function monentum(params::NamedTuple, coords::Matrix{Float64})
    mr = params.mr;
    mb = params.mb;
    l = params.l;
    mt = mr+3*mb;
    ptheta = mt * l^2*coords[:,2]/3;
    return ptheta;
end

t, coords = RKF45(sinPen, params, t0, tf, conds, epsilon, dtInit);
using PyPlot;
if (@isdefined(p1))
    PyPlot.close(p1)
end
p1 = PyPlot.figure(1)
PyPlot.plot(t, coords[:,1])
PyPlot.xlabel(L"t", rotation=0, fontsize=14)
PyPlot.ylabel(L"\theta",rotation=0,fontsize=14)
PyPlot.title("Single pendulum: angle with positive x-axis against time.", fontsize=16)
PyPlot.savefig("graphics/Figure 1 single pendulum theta vs t plot.png")
if (@isdefined(p2))
    PyPlot.close(p2)
end
p2 = PyPlot.figure(2)
PyPlot.plot(t, coords[:,2])
PyPlot.xlabel(L"t", rotation=0, fontsize=14)
PyPlot.ylabel(L"\dot{\theta}",rotation=0,fontsize=14)
PyPlot.title("Single pendulum: time derivative of angle with positive x-axis against time.", fontsize=16)
PyPlot.savefig("graphics/Figure 2 single pendulum theta dot vs t plot.png")
if (@isdefined(p3))
    PyPlot.close(p3)
end
p3 = PyPlot.figure(3)
PyPlot.plot(coords[:,1], coords[:,2])
PyPlot.xlabel(L"\theta", rotation=0, fontsize=14)
PyPlot.ylabel(L"\dot{\theta}", rotation=0, fontsize=14)
PyPlot.title("Single pendulum: phase space plot.", fontsize=16)
PyPlot.savefig("graphics/Figure 3 single pendulum theta dot vs theta plot.png")
x1 = params.l*cos.(coords[:,1]);
y1 = params.l*sin.(coords[:,1]);
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
        title = "Single pendulum",
        xlabel = "x", ylabel = "y",
        aspect = 1.0
    )

    # Initialise objects as vectors of 2D points
    line_data = [Point2f(0.0, 0.0), Point2f(x1_uni[1], y1_uni[1])]
    bob_data = [Point2f(x1_uni[1], y1_uni[1])]

    pendulum_line = lines!(ax, line_data, linewidth=2)
    bob = scatter!(ax, bob_data, color=:red, markersize=10)

    # Output path
    output_path = "graphics/single_pendulum.mp4"
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

animate(t, x1, y1, round(Int, tf*100))