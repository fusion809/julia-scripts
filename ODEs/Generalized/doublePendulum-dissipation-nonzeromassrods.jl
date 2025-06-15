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
m1r =  1.0, m1b =  1.0, m2r =  1.0, m2b =  1.0, 
b1b = 0.15, b1r = 0.15, b2b = 0.15, b2r = 0.15, 
c1b = 0.06, c1r = 0.06, c2b = 0.06, c2r = 0.06);

# Initial conditions and domain of integration
t0 = 0.0;
tf = 60.0;
theta10, dtheta10, theta20, dtheta20 = 0.0, 0.0, 0.0, 0.0;
conds = @SVector [theta10, dtheta10, theta20, dtheta20];

# Error tolerance and initial step size
epsilon = 1e-10;
dtInitial = 0.1;

# Solve problem and extract solution values
@time begin
t, vars = RKF45(doublePen, params, t0, tf, conds, epsilon, dtInitial);
end

# Extract dependent variables
theta1, dtheta1, theta2, dtheta2 = vars[:,1], vars[:,2], vars[:,3], vars[:,4];

# Positions of the pendulum bobs
x1 = params.r1*cos.(theta1);
y1 = params.r1*sin.(theta1);
x2 = x1 + params.r2*cos.(theta2);
y2 = y1 + params.r2*sin.(theta2);

# Plots
PyPlot.figure(1);
PyPlot.clf();
PyPlot.plot(t, theta1, label=L"\theta_1");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\theta_1");
PyPlot.legend();
PyPlot.savefig("graphics/Figure 1 Time plot of theta1.png")
# Time plot of dtheta1
PyPlot.figure(2);
PyPlot.clf();
PyPlot.plot(t, dtheta1, label=L"\dfrac{d\theta_1}{dt}");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\dfrac{d\theta_1}{dt}");
PyPlot.legend();
PyPlot.savefig("graphics/Figure 2 Time plot of dtheta1.png")
# Time plot of theta2
PyPlot.figure(3);
PyPlot.clf();
PyPlot.plot(t, theta2, label=L"\theta_2");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\theta_2");
PyPlot.legend();
PyPlot.savefig("graphics/Figure 3 Time plot of theta2.png")
# Time plot of dtheta2
PyPlot.figure(4);
PyPlot.clf();
PyPlot.plot(t, dtheta2, label=L"\dfrac{d\theta_2}{dt}");
PyPlot.xlabel(L"t");
PyPlot.ylabel(L"\dfrac{d\theta_2}{dt}");
PyPlot.legend();
PyPlot.savefig("graphics/Figure 4 Time plot of dtheta2.png")
# Phase plot theta2 vs theta1
PyPlot.figure(5);
PyPlot.clf();
PyPlot.plot(theta1, theta2, label=L"\theta_2");
PyPlot.xlabel(L"\theta_1");
PyPlot.ylabel(L"\theta_2");
PyPlot.legend();
PyPlot.savefig("graphics/Figure 5 Phase plot of theta2 against theta1.png")
# Phase plot dtheta1 vs theta1
PyPlot.figure(6);
PyPlot.clf();
PyPlot.plot(theta1, dtheta1, label=L"\dfrac{d\theta_1}{dt}");
PyPlot.xlabel(L"\theta_1");
PyPlot.ylabel(L"\dfrac{d\theta_1}{dt}");
PyPlot.legend();
PyPlot.savefig("graphics/Figure 6 Phase plot of dtheta1 against theta1.png")
# Phase plot dtheta2 vs theta2
PyPlot.figure(7);
PyPlot.clf();
PyPlot.plot(theta2, dtheta2, label=L"\dfrac{d\theta_2}{dt}");
PyPlot.xlabel(L"\theta_2");
PyPlot.ylabel(L"\dfrac{d\theta_2}{dt}");
PyPlot.legend();
PyPlot.savefig("graphics/Figure 7 Phase plot of dtheta2 against theta2.png")
# All angles and derivatives
PyPlot.figure(8);
PyPlot.clf();
PyPlot.plot(t, vars[:,1], label=L"\theta_1");
PyPlot.plot(t, vars[:,2], label=L"\dfrac{d\theta_1}{dt}");
PyPlot.plot(t, vars[:,3], label=L"\theta_2");
PyPlot.plot(t, vars[:,4], label=L"\dfrac{d\theta_2}{dt}");
PyPlot.xlabel(L"t")
PyPlot.legend()
PyPlot.savefig("graphics/Figure 8 All angles and derivatives against time.png")
# Bob #1
PyPlot.figure(9);
PyPlot.clf();
PyPlot.plot(x1, y1, label="Pendulum bob 1 location")
PyPlot.legend()
PyPlot.savefig("graphics/Figure 9 Pendulum bob 1 location.png")
# Bob #2
PyPlot.figure(10);
PyPlot.clf();
PyPlot.plot(x2, y2, label="Pendulum bob 2 location")
PyPlot.legend()
PyPlot.savefig("graphics/Figure 10 Pendulum bob 2 location.png")
# Bob #1 and #2
PyPlot.figure(11);
PyPlot.clf();
PyPlot.plot(x1, y1, label="Pendulum bob 1 location")
PyPlot.plot(x2, y2, label="Pendulum bob 2 location")
PyPlot.legend()
PyPlot.savefig("graphics/Figure 11 Pendulum bob 1 and 2 location.png")

# Setup figure and scene
using CairoMakie;
line1_data = Observable([[0.0, x1[1]], [0.0, y1[1]]])
line2_data = Observable([[x1[1], x2[1]], [y1[1], y2[1]]])
mass_data  = Observable((x=[x1[1], x2[1]], y=[y1[1], y2[1]]))
time_text  = Observable("t = 0.00 s")
f = CairoMakie.Figure(resolution = (600, 600))
ax = Axis(f[1, 1], aspect = 1, limits = (-2.2*(params.r1+params.r2)/2, 2.2*(params.r1+params.r2)/2, -2.2*(params.r1+params.r2)/2, 2.2*(params.r1+params.r2)/2))
pendulum_line1 = lines!(ax, line1_data[][1], line1_data[][2], color = :blue)
pendulum_line2 = lines!(ax, line2_data[][1], line2_data[][2], color = :red)
masses = scatter!(ax, mass_data[].x, mass_data[].y, color = [:blue, :red], markersize = 10)
t_uni = LinRange(t0, tf, 10001);
splx1 = Spline1D(t, x1)
splx2 = Spline1D(t, x2)
sply1 = Spline1D(t, y1)
sply2 = Spline1D(t, y2)
x1_uni = evaluate(splx1, t_uni)
x2_uni = evaluate(splx2, t_uni)
y1_uni = evaluate(sply1, t_uni)
y2_uni = evaluate(sply2, t_uni)
text!(ax, time_text, position = Point2f(0, 2.1), align = (:left, :top), fontsize = 18, color = :black)

# Create and save animation
Dt = step(t_uni)
record(f, "graphics/Figure 12 Double pendulum animation.mp4", 1:length(t_uni); framerate = round(Int, 1/Dt)) do i
    pendulum_line1[1] = [0, x1_uni[i]]
    pendulum_line1[2] = [0, y1_uni[i]]
    pendulum_line2[1] = [x1_uni[i], x2_uni[i]]
    pendulum_line2[2] = [y1_uni[i], y2_uni[i]]
    masses[1] = [x1_uni[i], x2_uni[i]]
    masses[2] = [y1_uni[i], y2_uni[i]]
    time_text[]  = "t = $(round(t_uni[i], digits = 2)) s"
end
