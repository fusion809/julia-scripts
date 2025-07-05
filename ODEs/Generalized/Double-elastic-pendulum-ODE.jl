include("RKF45.jl")
using LinearAlgebra;
function elastic_double_pendulum!(params::NamedTuple, t::Float64, u::SVector{8,Float64})::SVector{8,Float64}
    # Unpack variables
    m1 = params.m1;
    m2 = params.m2;
    l1 = params.l1;
    l2 = params.l2;
    k1 = params.k1;
    k2 = params.k2;
    b1 = params.b1;
    c1 = params.c1;
    b2 = params.b2;
    c2 = params.c2;
    g  = params.g;
    r1   = u[1];
    dr1  = u[2];
    r2   = u[3];
    dr2  = u[4];
    th1  = u[5];
    dth1 = u[6];
    th2  = u[7];
    dth2 = u[8];
    del = th2 - th1
    v1  = sqrt(dr1^2+r1^2*dth1^2);
    v2  = sqrt(v1^2 + dr2^2 + r2^2*dth2^2 + 2*cos(del)*(dr1*dr2 + r1*r2*dth1*dth2)+2*sin(del)*(r1*dr2*dth1-dr1*r2*dth2));
    Qr1 = -(b1 + c1*v1)*dr1 - (b2 + c2*v2)*(dr1 + dr2*cos(del)-r2*dth2*sin(del));
    Qr2 = -(b2 + c2*v2)*(dr1*cos(del) + r1*dth1*sin(del)+dr2);
    Qt1 = -(b1 + c1*v1)*r1^2*dth1 - (b2 + c2*v2)*(r1^2*dth1 + r1*dr2*sin(del) + r1*r2*dth2*cos(del))
    Qt2 = -(b2 + c2*v2)*(r2^2*dth2 - dr1*r2*sin(del) + r1*r2*dth1*cos(del))
    # Force equations (from your full system)
    # Use the original non-rearranged equations here, e.g.:

    A = Matrix{Float64}(I, 4, 4);
    b = zeros(Float64, 4,1);
    A[1,2] = m2/((m1+m2))*cos(del);
    A[1,4] = -m2*r2/(m1+m2)*sin(del);
    A[2,1] = cos(del);
    A[2,3] = r1*sin(del);
    A[3,2] = m2/((m1+m2)*r1)*sin(del);
    A[3,4] = m2*r2/((m1+m2)*r1)*cos(del);
    A[4,1] = -sin(del)/r2;
    A[4,3] = cos(del)*r1/r2;
    b[1] = r1*dth1^2 - g*sin(th1) + m2/(m1+m2)*(cos(del)*r2*dth2^2+sin(del)*2*dr2*dth2) + (Qr1 - k1*(r1-l1))/(m1+m2);
    b[2] = r2*dth2^2 - g*sin(th2) + cos(del)*r1*dth1^2-2*sin(del)*dr1*dth1 + (Qr2-k2*(r2-l2))/m2;
    b[3] = -2*dr1*dth1/r1 - g*cos(th1)/r1 - m2/((m1+m2)*r1)*(cos(del)*2*dr2*dth2-sin(del)*r2*dth2^2) + Qt1/((m1+m2)*r1^2);
    b[4] = -2*dr2*dth2/r2 - g*cos(th2)/r2 - cos(del)/r2*2*dr1*dth1 - sin(del)*r1*dth1^2/r2 + Qt2/(m2*r2^2);
    d2m = A\b;
    return [dr1, d2m[1], dr2, d2m[2], dth1, d2m[3], dth2, d2m[4]];
end

# Initial conditions
u0 = @SVector [1.0, 0.0, 1, 0, 0, 0, 0, 0]

params = (m1=1.0, m2=1.0, g=9.81, k1=1.0, k2=1.0, l1=1.0, l2=1.0, b1=0.0, c1=0.0, b2=0.0, c2=0.0)

#tspan = (0.0, 1.0)
t0 = 0.0;
tf = 300.0;
epsilon = 1e-9;
dtInitial = 1e-3;
t, vars = RKF45(elastic_double_pendulum!, params, t0, tf, u0, epsilon, dtInitial)
r1   = vars[:,1];
dr1  = vars[:,2];
r2   = vars[:,3];
dr2  = vars[:,4];
th1  = vars[:,5];
dth1 = vars[:,6];
th2  = vars[:,7];
dth2 = vars[:,8];
using PyPlot;
if @isdefined(p1)
        PyPlot.close(p1)
end
p1 = PyPlot.figure(1)
PyPlot.plot(t, r1)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"r_1", rotation=0, fontsize=14)
PyPlot.title("Figure 1 Double elastic pendulum \$r_1\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Figure 1 Double elastic pendulum r1 vs t.png")
if @isdefined(p2)
        PyPlot.close(p2)
end
p2 = PyPlot.figure(2)
PyPlot.plot(t, r2)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"r_2", rotation=0, fontsize=14)
PyPlot.title("Figure 2 Double elastic pendulum \$r_2\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Figure 2 Double elastic pendulum r2 vs t.png")
if @isdefined(p3)
        PyPlot.close(p3)
end
p3 = PyPlot.figure(3)
PyPlot.plot(t, dr1)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\dot{r}_1", rotation=0, fontsize=14)
PyPlot.title("Figure 3 Double elastic pendulum \$\\dot{r}_1\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Figure 3 Double elastic pendulum r1dot vs t.png")
if @isdefined(p4)
        PyPlot.close(p4)
end
p4 = PyPlot.figure(4)
PyPlot.plot(t, dr2)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\dot{r}_2", rotation=0, fontsize=14)
PyPlot.title("Figure 4 Double elastic pendulum \$\\dot{r}_2\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Figure 4 Double elastic pendulum r2dot vs t.png")
if @isdefined(p5)
        PyPlot.close(p5)
end
p5 = PyPlot.figure(5)
PyPlot.plot(r1, dr1)
PyPlot.xlabel(L"r_1", fontsize=14)
PyPlot.ylabel(L"\dot{r}_1", rotation=0, fontsize=14)
PyPlot.title("Figure 5 Double elastic pendulum \$\\dot{r}_1\$ vs \$r_1\$.", fontsize=16)
PyPlot.savefig("graphics/Figure 5 Double elastic pendulum r1dot vs r1.png")
if @isdefined(p6)
        PyPlot.close(p6)
end
p6 = PyPlot.figure(6)
PyPlot.plot(r2, dr2)
PyPlot.xlabel(L"r_2", fontsize=14)
PyPlot.ylabel(L"\dot{r}_2", rotation=0, fontsize=14)
PyPlot.title("Figure 6 Double elastic pendulum \$\\dot{r}_2\$ vs \$r_2\$.", fontsize=16)
PyPlot.savefig("graphics/Figure 6 Double elastic pendulum r2dot vs r2.png")
if @isdefined(p7)
        PyPlot.close(p7)
end
p7 = PyPlot.figure(7)
PyPlot.plot(t, th1)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\theta_1", rotation=0, fontsize=14)
PyPlot.title("Figure 7 Double elastic pendulum \$\\theta_1\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Figure 7 Double elastic pendulum theta1 vs t.png")
if @isdefined(p8)
        PyPlot.close(p8)
end
p8 = PyPlot.figure(8)
PyPlot.plot(t, th2)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\theta_2", rotation=0, fontsize=14)
PyPlot.title("Figure 8 Double elastic pendulum \$\\theta_2\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Figure 8 Double elastic pendulum theta2 vs t.png")
if @isdefined(p9)
        PyPlot.close(p9)
end
p9 = PyPlot.figure(9)
PyPlot.plot(t, dth1)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}_1", rotation=0, fontsize=14)
PyPlot.title("Figure 9 Double elastic pendulum \$\\dot{\\theta}_1\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Figure 9 Double elastic pendulum theta1dot vs t.png")
if @isdefined(p10)
        PyPlot.close(p10)
end
p10 = PyPlot.figure(10)
PyPlot.plot(t, dth2)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}_2", rotation=0, fontsize=14)
PyPlot.title("Figure 10 Double elastic pendulum \$\\dot{\\theta}_2\$ vs t.", fontsize=16)
PyPlot.savefig("graphics/Figure 10 Double elastic pendulum theta2dot vs t.png")
if @isdefined(p11)
        PyPlot.close(p11)
end
p11 = PyPlot.figure(11)
PyPlot.plot(th1, dth1)
PyPlot.xlabel(L"\theta_1", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}_1", rotation=0, fontsize=14)
PyPlot.title("Figure 11 Double elastic pendulum \$\\dot{\\theta}_1\$ vs \$\\theta_1\$.", fontsize=16)
PyPlot.savefig("graphics/Figure 11 Double elastic pendulum theta1dot vs theta1.png")
if @isdefined(p12)
        PyPlot.close(p12)
end
p12 = PyPlot.figure(12)
PyPlot.plot(th2, dth2)
PyPlot.xlabel(L"\theta_2", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}_2", rotation=0, fontsize=14)
PyPlot.title("Figure 12 Double elastic pendulum \$\\dot{\\theta}_2\$ vs \$\\theta_2\$.", fontsize=16)
PyPlot.savefig("graphics/Figure 12 Double elastic pendulum theta2dot vs theta2.png")

using CairoMakie, Dierckx;
function pendpos(r1, r2, theta1, theta2)
        if (length(r1) == 1)
                x1 = r1*cos.(th1);
                y1 = r1*sin.(th1);
        else
                x1 = r1.*cos.(th1);
                y1 = r1.*sin.(th1);
        end
        if (length(r2) == 1)
                x2 = x1 .+ r2*cos.(th2);
                y2 = y1 .+ r2*sin.(th2);
        else
                x2 = x1 .+ r2.*cos.(th2);
                y2 = y1 .+ r2.*sin.(th2);
        end
        return x1, y1, x2, y2;
end

function spline_double_pendulum(t, t_uni, x1, y1, x2, y2)
        splx1          = Spline1D(t, x1)
        splx2          = Spline1D(t, x2)
        sply1          = Spline1D(t, y1)
        sply2          = Spline1D(t, y2)
        x1_uni         = evaluate(splx1, t_uni)
        x2_uni         = evaluate(splx2, t_uni)
        y1_uni         = evaluate(sply1, t_uni)
        y2_uni         = evaluate(sply2, t_uni)
        return x1_uni, y1_uni, x2_uni, y2_uni
end

function double_pendulum_anim(t, r1, r2, theta1, theta2, N=10001, padrat = 0.1, filename="Figure 13 Double elastic pendulum animation.mp4")
        x1, y1, x2, y2 = pendpos(r1, r2, theta1, theta2);
        # Setup animation
        tf             = maximum(t)
        t0             = minimum(t)
        dt             = (tf-t0)/N;
        t_uni          = t0:dt:tf;
        if (length(r1) > 1)
                padding = padrat*maximum(r1 .+ r2);
        else
                padding = padrat*(r1+r2);
        end
        line1_data     = Observable([[0.0, x1[1]], [0.0, y1[1]]])
        line2_data     = Observable([[x1[1], x2[1]], [y1[1], y2[1]]])
        mass_data      = Observable((x=[x1[1], x2[1]], y=[y1[1], y2[1]]))
        time_text      = Observable("t = 0.00 s")
        f              = CairoMakie.Figure(resolution = (600, 600))
        ax             = Axis(f[1, 1], aspect = 1, 
        limits = (minimum([x1 x2])-padding, maximum([x1 x2])+padding, 
        minimum([y1 y2])-padding, maximum([y1 y2])+padding))
        pendulum_line1 = lines!(ax, line1_data[][1], line1_data[][2], color = :blue)
        pendulum_line2 = lines!(ax, line2_data[][1], line2_data[][2], color = :red)
        masses         = scatter!(ax, mass_data[].x, mass_data[].y, color = [:blue, :red], markersize = 10)
        text!(ax, time_text, position = Point2f(0, 2.1), align = (:left, :top), 
        fontsize = 18, color = :black)
        x1_uni, y1_uni, x2_uni, y2_uni = spline_double_pendulum(t, t_uni, x1, y1, x2, y2)

        # Create and save animation
        record(f, "graphics/$filename", 1:N; framerate = round(Int, 1/dt)) do i
        pendulum_line1[1] = [0, x1_uni[i]]
        pendulum_line1[2] = [0, y1_uni[i]]
        pendulum_line2[1] = [x1_uni[i], x2_uni[i]]
        pendulum_line2[2] = [y1_uni[i], y2_uni[i]]
        masses[1]         = [x1_uni[i], x2_uni[i]]
        masses[2]         = [y1_uni[i], y2_uni[i]]
        time_text[]       = "t = $(round(t_uni[i], digits = 2)) s"
        end
end
double_pendulum_anim(t, r1, r2, th1, th2)