using PyPlot, FFMPEG, StaticArrays, Dierckx;
include("RKF45.jl");

function doublePen(params::NamedTuple, t::Float64, vars::SVector{4, Float64})::SVector{4, Float64}
    # Extract parameters
    g   = params.g;
    r1  = params.r1;
    r2  = params.r2;
    m1b = params.m1b;
    m1r = params.m1r;
    m2b = params.m2b;
    m2r = params.m2r;
    b1b = params.b1b;
    b1r = params.b1r;
    b2b = params.b2b;
    b2r = params.b2r;
    c1b = params.c1b;
    c1r = params.c1r;
    c2b = params.c2b;
    c2r = params.c2r;

    # Dependent variables
    theta1   = vars[1];
    dtheta1  = vars[2];
    theta2   = vars[3];
    dtheta2  = vars[4];

    # Variables to simplify later definitions
    outercoef = 1/((m2r/12+m2b)*r2 - (m2b^2*r2*cos(theta1-theta2)^2)/(m1r/12 + m1b+m2b));
    mass = m1r/2 + m2r + m1b + m2b; 
    v1b = r1*dtheta1;
    v1r = v1b/2;
    v2b = sqrt(r1^2*dtheta1^2+r2^2*dtheta2^2+2*r1*r2*dtheta1*dtheta2*cos(theta1-theta2));
    v2r = sqrt(r1^2*dtheta1^2+(r2^2*dtheta2^2)/4 + r2*dtheta1*dtheta2*cos(theta1-theta2));
    drag1b = (b1b + c1b*v1b)*v1b;
    drag1r = (b1r + c1r*v1r)*v1r/2;
    drag2b = (b2b + c2b*v2b)*(r1*dtheta1 + r2*dtheta2*cos(theta1-theta2));
    drag2b2 = (b2b + c2b*v2b)*(r1*dtheta1*cos(theta1-theta2) + r2*dtheta2);
    drag2r = (b2r + c2r*v2r)*(r1*dtheta1 + r2*dtheta2*cos(theta1-theta2)/2);
    drag2r2 = 1/4*(b2r + c2r*v2r)*(2*r1*dtheta1*cos(theta1-theta2) + r2*dtheta2);
    innercoef = -m2b*cos(theta1-theta2)/(m1r/12+m1b+m2b);
    inner = innercoef*(-m2b*r2*dtheta2^2*sin(theta1-theta2) - g*cos(theta1)*(m1r/2+m2r+m1b+m2b)-drag1b-drag2b-drag1r-drag2r);
    extra = m2b*(r1*dtheta1^2*sin(theta1-theta2)-g*cos(theta2))-drag2r2-drag2b2;
    d2theta2 = outercoef*(inner+extra);
    outercoef1 = 1/((m1r/12+m1b+m2b)*r1);
    inner11 = -m2b*r2*(d2theta2*cos(theta1-theta2)+dtheta2^2*sin(theta1-theta2));
    inner12 = -g*cos(theta1)*mass - drag1b - drag2b - drag1r - drag2r;
    d2theta1 = outercoef1*(inner11 + inner12);

    # Return statement
    return [dtheta1, d2theta1, dtheta2, d2theta2];
end

# Problem parameters
params = (g = 9.81, r1 = 1.0, r2 = 1.0, 
m1r =  0.0, m1b =  1.0, m2r =  0.0, m2b =  1.0, 
b1b = 0.00, b1r = 0.00, b2b = 0.00, b2r = 0.00, 
c1b = 0.00, c1r = 0.00, c2b = 0.00, c2r = 0.00);
# params = (g = 9.81, r1 = 1.0, r2 = 1.0, 
# m1r =  1.0, m1b =  1.0, m2r =  1.0, m2b =  1.0, 
# b1b = 0.15, b1r = 0.15, b2b = 0.15, b2r = 0.15, 
# c1b = 0.06, c1r = 0.06, c2b = 0.06, c2r = 0.06);

# Initial conditions and domain of integration
t0 = 0.0;
tf = 60.0;
theta10, dtheta10, theta20, dtheta20 = 0.0, 0.0, 0.0, 0.0;
theta11, dtheta11, theta21, dtheta21 = 0.001*pi/180, 0.0, 0.0, 0.0;
conds0 = @SVector [theta10, dtheta10, theta20, dtheta20];
conds1 = @SVector [theta11, dtheta11, theta21, dtheta21];

# Error tolerance and initial step size
epsilon = 1e-9;
dtInitial = 0.1;

# Solve problem and extract solution values
@time begin
t_0, vars0 = RKF45(doublePen, params, t0, tf, conds0, epsilon, dtInitial);
end
@time begin
t_1, vars1 = RKF45(doublePen, params, t0, tf, conds1, epsilon, dtInitial);
end

# Extract dependent variables
theta1_0, dtheta1_0, theta2_0, dtheta2_0 = vars0[:,1], vars0[:,2], vars0[:,3], vars0[:,4];
theta1_1, dtheta1_1, theta2_1, dtheta2_1 = vars1[:,1], vars1[:,2], vars1[:,3], vars1[:,4];

# Positions of the pendulum bobs
x1_0 = params.r1*cos.(theta1_0);
y1_0 = params.r1*sin.(theta1_0);
x2_0 = x1_0 + params.r2*cos.(theta2_0);
y2_0 = y1_0 + params.r2*sin.(theta2_0);
x1_1 = params.r1*cos.(theta1_1);
y1_1 = params.r1*sin.(theta1_1);
x2_1 = x1_1 + params.r2*cos.(theta2_1);
y2_1 = y1_1 + params.r2*sin.(theta2_1);

# Plots
PyPlot.figure(1);
PyPlot.clf();
PyPlot.plot(t_0, theta1_0, label=L"\theta_{1,0}");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\theta_{1,0}");
PyPlot.legend();
PyPlot.title("\$\\theta_{1,0}\$ vs \$t\$")
PyPlot.savefig("graphics/Figure 1 Time plot of theta1_0.png")
# Time plot of dtheta1
PyPlot.figure(2);
PyPlot.clf();
PyPlot.plot(t_0, dtheta1_0, label=L"\dot{\theta}_{1,0}");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\dot{\theta}_{1,0}");
PyPlot.legend();
PyPlot.title("\$\\dot{\\theta}_{1,0}\$ vs \$t\$")
PyPlot.savefig("graphics/Figure 2 Time plot of dtheta1_0.png")
# Time plot of theta2
PyPlot.figure(3);
PyPlot.clf();
PyPlot.plot(t_0, theta2_0, label=L"\theta_{2,0}");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\theta_{2,0}");
PyPlot.legend();
PyPlot.title("\$\\theta_{2,0}\$ vs \$t\$")
PyPlot.savefig("graphics/Figure 3 Time plot of theta2_0.png")
# Time plot of dtheta2
PyPlot.figure(4);
PyPlot.clf();
PyPlot.plot(t_0, dtheta2_0, label=L"\dot{\theta}_{2,0}");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\dot{\theta}_{2,0}");
PyPlot.legend();
PyPlot.title("\$\\dot{\\theta}_{2,0}\$ vs \$t\$")
PyPlot.savefig("graphics/Figure 4 Time plot of dtheta2_0.png")
# Phase plot theta2 vs theta1
PyPlot.figure(5);
PyPlot.clf();
PyPlot.plot(theta1_0, theta2_0, label=L"\theta_{2,0}");
PyPlot.xlabel(L"\theta_{1,0}");
PyPlot.ylabel(L"\theta_{2,0}");
PyPlot.legend();
PyPlot.title("\$\\theta_{2,0}\$ vs \$\\theta_{1,0}\$")
PyPlot.savefig("graphics/Figure 5 Phase plot of theta2_0 against theta1_0.png")
# Phase plot dtheta1 vs theta1
PyPlot.figure(6);
PyPlot.clf();
PyPlot.plot(theta1_0, dtheta1_0, label=L"\dot{\theta}_{1,0}");
PyPlot.xlabel(L"\theta_{1,0}");
PyPlot.ylabel(L"\dot{\theta}_{1,0}");
PyPlot.legend();
PyPlot.title("\$\\dot{\\theta}_{1,0}\$ vs \$\\theta_{1,0}\$")
PyPlot.savefig("graphics/Figure 6 Phase plot of dtheta1_0 against theta1_0.png")
# Phase plot dtheta2 vs theta2
PyPlot.figure(7);
PyPlot.clf();
PyPlot.plot(theta2_0, dtheta2_0, label=L"\dot{\theta}_{2,0}");
PyPlot.xlabel(L"\theta_{2,0}");
PyPlot.ylabel(L"\dot{\theta}_{2,0}");
PyPlot.legend();
PyPlot.title("\$\\dot{\\theta}_{2,0}\$ vs \$\\theta_{2,0}\$")
PyPlot.savefig("graphics/Figure 7 Phase plot of dtheta2_0 against theta2_0.png")
# All angles and derivatives
PyPlot.figure(8);
PyPlot.clf();
PyPlot.plot(t_0, vars0[:,1], label=L"\theta_{1,0}");
PyPlot.plot(t_0, vars0[:,2], label=L"\dot{\theta}_{1,0}");
PyPlot.plot(t_0, vars0[:,3], label=L"\theta_{2,0}");
PyPlot.plot(t_0, vars0[:,4], label=L"\dot{\theta}_{2,0}");
PyPlot.xlabel(L"t")
PyPlot.legend()
PyPlot.title("All angles and derivatives for first init conds against time")
PyPlot.savefig("graphics/Figure 8 All angles and derivatives for first init conds against time.png")
# Time plot of theta1
PyPlot.figure(9);
PyPlot.clf();
PyPlot.plot(t_1, theta1_1, label=L"\theta_{1,1}");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\theta_{1,1}");
PyPlot.legend();
PyPlot.title("\$\\theta_{1,1}\$ vs \$t\$")
PyPlot.savefig("graphics/Figure 9 Time plot of theta1_1.png")
# Time plot of dtheta1
PyPlot.figure(10);
PyPlot.clf();
PyPlot.plot(t_1, dtheta1_1, label=L"\dot{\theta}_{1,1}");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\dot{\theta}_{1,1}");
PyPlot.legend();
PyPlot.title("\$\\dot{\\theta}_{1,1}\$ vs \$t\$")
PyPlot.savefig("graphics/Figure 10 Time plot of dtheta1_1.png")
# Time plot of theta2
PyPlot.figure(11);
PyPlot.clf();
PyPlot.plot(t_1, theta2_1, label=L"\theta_{2,1}");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\theta_{2,1}");
PyPlot.legend();
PyPlot.title("\$\\theta_{2,1}\$ vs \$t\$")
PyPlot.savefig("graphics/Figure 11 Time plot of theta2_1.png")
# Time plot of dtheta2
PyPlot.figure(12);
PyPlot.clf();
PyPlot.plot(t_1, dtheta2_1, label=L"\dot{\theta}_{2,1}");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\dot{\theta}_{2,1}");
PyPlot.legend();
PyPlot.title("\$\\dot{\\theta}_{2,1}\$ vs \$t\$")
PyPlot.savefig("graphics/Figure 12 Time plot of dtheta2_1.png")
# Phase plot theta2 vs theta1
PyPlot.figure(13);
PyPlot.clf();
PyPlot.plot(theta1_1, theta2_1, label=L"\theta_{2,1}");
PyPlot.xlabel(L"\theta_{1,1}");
PyPlot.ylabel(L"\theta_{2,1}");
PyPlot.legend();
PyPlot.title("\$\\theta_{2,1}\$ vs \$\\theta_{1,1}\$")
PyPlot.savefig("graphics/Figure 13 Phase plot of theta2_1 against theta1_1.png")
# Phase plot dtheta1 vs theta1
PyPlot.figure(14);
PyPlot.clf();
PyPlot.plot(theta1_1, dtheta1_1, label=L"\dot{\theta}_{1,1}");
PyPlot.xlabel(L"\theta_{1,1}");
PyPlot.ylabel(L"\dot{\theta}_{1,1}");
PyPlot.legend();
PyPlot.title("\$\\dot{\\theta}_{1,1}\$ vs \$\\theta_{1,1}\$")
PyPlot.savefig("graphics/Figure 14 Phase plot of dtheta1_1 against theta1_1.png")
# Phase plot dtheta2 vs theta2
PyPlot.figure(15);
PyPlot.clf();
PyPlot.plot(theta2_1, dtheta2_1, label=L"\dot{\theta}_{2,1}");
PyPlot.xlabel(L"\theta_{2,1}");
PyPlot.ylabel(L"\dot{\theta}_{2,1}");
PyPlot.legend();
PyPlot.title("\$\\dot{\\theta}_{2,1}\$ vs \$\\theta_{2,1}\$")
PyPlot.savefig("graphics/Figure 15 Phase plot of dtheta2_1 against theta2_1.png")
# All angles and derivatives
PyPlot.figure(16);
PyPlot.clf();
PyPlot.plot(t_1, vars1[:,1], label=L"\theta_{1,1}");
PyPlot.plot(t_1, vars1[:,2], label=L"\dot{\theta}_{1,1}");
PyPlot.plot(t_1, vars1[:,3], label=L"\theta_{2,1}");
PyPlot.plot(t_1, vars1[:,4], label=L"\dot{d\theta}_{2,1}");
PyPlot.xlabel(L"t")
PyPlot.title("All angles and derivatives for second init conds against time")
PyPlot.legend("Figure 16 All angles and derivatives for second init conds against time")
PyPlot.savefig("graphics/Figure 16 All angles and derivatives for second init conds against time.png")

# Setup animation
t_uni = LinRange(t0, tf, 6001);
splx1_0 = Spline1D(t_0, x1_0)
splx2_0 = Spline1D(t_0, x2_0)
sply1_0 = Spline1D(t_0, y1_0)
sply2_0 = Spline1D(t_0, y2_0)
splx1_1 = Spline1D(t_1, x1_1)
splx2_1 = Spline1D(t_1, x2_1)
sply1_1 = Spline1D(t_1, y1_1)
sply2_1 = Spline1D(t_1, y2_1)
x1_uni_0 = evaluate(splx1_0, t_uni)
x2_uni_0 = evaluate(splx2_0, t_uni)
y1_uni_0 = evaluate(sply1_0, t_uni)
y2_uni_0 = evaluate(sply2_0, t_uni)
x1_uni_1 = evaluate(splx1_1, t_uni)
x2_uni_1 = evaluate(splx2_1, t_uni)
y1_uni_1 = evaluate(sply1_1, t_uni)
y2_uni_1 = evaluate(sply2_1, t_uni)

using CairoMakie

CairoMakie.activate!()

# Time steps
N = length(t_uni)

# Set up figure and axis
f = CairoMakie.Figure(resolution = (600, 600))
ax = Axis(f[1, 1], xlabel = "x", ylabel = "y", aspect = DataAspect())
ax.title = "Double pendulums simulation, θ₁(t=0) differ by $(round(180/pi*(theta11-theta10), digits=4)) degrees."
ax.xlabel = L"x"
ax.ylabel = L"y"

# Frame 1 initial positions
pts_line0 = [Point2f(0, 0), Point2f(x1_0[1], y1_0[1]), Point2f(x2_0[1], y2_0[1])]
pts_line1 = [Point2f(0, 0), Point2f(x1_1[1], y1_1[1]), Point2f(x2_1[1], y2_1[1])]

# Pendulum 0 (blue)
pendulum_line0 = lines!(ax, pts_line0, color = :blue)
bob1_0 = scatter!(ax, [Point2f(x1_0[1], y1_0[1])], markersize = 10, color = :blue)
bob2_0 = scatter!(ax, [Point2f(x2_0[1], y2_0[1])], markersize = 10, color = :blue)

# Pendulum 1 (red)
pendulum_line1 = lines!(ax, pts_line1, color = :red)
bob1_1 = scatter!(ax, [Point2f(x1_1[1], y1_1[1])], markersize = 10, color = :red)
bob2_1 = scatter!(ax, [Point2f(x2_1[1], y2_1[1])], markersize = 10, color = :red)

# Set axis limits
x_all = [x1_0; x2_0; x1_1; x2_1]
y_all = [y1_0; y2_0; y1_1; y2_1]
xlims!(ax, minimum(x_all) - 0.5, maximum(x_all) + 0.5)
ylims!(ax, minimum(y_all) - 0.5, maximum(y_all) + 0.5)

# Animate
Dt = step(t_uni)
record(f, "graphics/Figure 17 Double pendulums masslessrods no dissipation.mp4", 1:N; framerate = round(Int, 1/Dt)) do i
    pendulum_line0[1][] = [Point2f(0, 0), Point2f(x1_uni_0[i], y1_uni_0[i]), 
    Point2f(x2_uni_0[i], y2_uni_0[i])]
    bob1_0[1][] = [Point2f(x1_uni_0[i], y1_uni_0[i])]
    bob2_0[1][] = [Point2f(x2_uni_0[i], y2_uni_0[i])]

    pendulum_line1[1][] = [Point2f(0, 0), Point2f(x1_uni_1[i], y1_uni_1[i]), 
    Point2f(x2_uni_1[i], y2_uni_1[i])]
    bob1_1[1][] = [Point2f(x1_uni_1[i], y1_uni_1[i])]
    bob2_1[1][] = [Point2f(x2_uni_1[i], y2_uni_1[i])]
    title
end