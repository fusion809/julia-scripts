include("RKF45.jl")
using LinearAlgebra
params   = (m1 = 1, m2=1, g=9.81, r=1, k=1, l=1, 
b1=0.0, c1=0.00, b2=0.0, c2=0.00);
z0       = 0.0;
dz0      = 0.0;
theta10  = 0;
dtheta10 = 0;
theta20  = 0;
dtheta20 = 0;
t0       = 0.0;
tf       = 200.0;
conds    = @SVector [z0, dz0, theta10, dtheta10, theta20, dtheta20];
dtInit   = 1e-3;
epsilon  = 1e-9;

function DP2E(params::NamedTuple, t::Float64, coords::SVector{6, Float64})::SVector{6,Float64}
    g  = params.g;
    l  = params.l;
    r  = params.r;
    m1 = params.m1;
    m2 = params.m2;
    b1 = params.b1;
    b2 = params.b2;
    c1 = params.c1;
    c2 = params.c2;
    k  = params.k;
    theta1   = coords[3];
    dtheta1  = coords[4];
    theta2   = coords[5];
    dtheta2  = coords[6];
    z        = coords[1];
    dz       = coords[2];
    v1       = r*abs(dtheta1);
    v2       = sqrt(r^2*dtheta1^2 + dz^2 + (l+z)^2*dtheta2^2 + 2*r*dtheta1 * dz * sin(theta2-theta1) + 2*r*(l+z)*dtheta1*dtheta2*cos(theta2-theta1));
    bob1dr   = -(b1 + c1*v1);
    bob2dr   = -(b2 + c2*v2);
    Qz       = bob2dr * (r*dtheta1 * sin(theta2-theta1) + dz);
    Qtheta1  = bob1dr * (r^2*dtheta1) + bob2dr * (r^2*dtheta1 + r*dz*sin(theta2-theta1) + r*(l+z)*dtheta2*cos(theta2-theta1));
    Qtheta2  = bob2dr * (l+z) * (r*dtheta1*cos(theta2-theta1) + (l+z)*dtheta2);
    # firstout = -(m1+m2*(cos(theta2-theta1))^2)/m1;
    # d2theta2 = firstout * (Qtheta2/(m2*(l+z)^2)-g*cos(theta2)/(l+z)) - (r*dtheta1^2*sin(theta2-theta1)+2*dz*dtheta2)/(l+z) - (cos(theta2-theta1))/(m1*r*(l+z)) * (m2*g*r*sin(theta2)*sin(theta2-theta1) + Qtheta1 + r*(k*z-Qz)*sin(theta2-theta1)-(m1+m2)*g*r*cos(theta1));
    # outer    = 1/((m1+m2*(cos(theta2-theta1))^2)*r^2);
    # inner    = Qtheta1 - g*(m1+m2)*r*cos(theta1) + r*(k*z-Qz)*sin(theta2-theta1) - m2*r*((r*dtheta1^2*cos(theta2-theta1))*sin(theta2-theta1) + (2*dz*dtheta2 + (l+z)*dtheta2)*cos(theta2-theta1));
    # d2theta1 = outer * inner;
    # d2z      = -r*d2theta1 * sin(theta2-theta1) + r*dtheta1^2*cos(theta1-theta1) + (l+z)*dtheta2^2 - g*sin(theta2) + (Qz-k*z)/m2;
    A = Matrix{Float64}(I, 3, 3);
    b = zeros(3, 1);
    A[1, 2] = r*sin(theta2-theta1);
    b[1] = r*dtheta1^2*cos(theta2-theta1) + (l+z)*dtheta2^2 - g*sin(theta2) + (Qz-k*z)/m2;
    A[2, 1] = m2/((m1+m2)*r)*sin(theta2-theta1);
    A[2, 3] = m2/((m1+m2)*r)*(l+z)*cos(theta2-theta1);
    b[2] = -g/r*cos(theta1) + Qtheta2/((m1+m2)*r^2) - m2/((m1+m2)*r)*(-(l+z)*dtheta2^2*sin(theta2-theta1)+2*dz*dtheta2*cos(theta2-theta1));
    A[3, 2] = r*cos(theta2-theta1)/(l+z);
    b[3] = Qtheta2/(m2*(l+z)^2) - 1/(l+z)*(2*dz*dtheta2 + r*dtheta1^2*sin(theta2-theta1)+g*cos(theta2));
    secder = A\b;
    return [dz, secder[1], dtheta1, secder[2], dtheta2, secder[3]];
end

t, coords = RKF45(DP2E, params, t0, tf, conds, epsilon, dtInit, "relative", 1e-6);
theta1    = coords[:,3];
dtheta1   = coords[:,4];
theta2    = coords[:,5];
dtheta2   = coords[:,6];
z         = coords[:,1];
dz        = coords[:,2];
using PyPlot
if @isdefined(p1)
    PyPlot.close(p1)
end
p1 = PyPlot.figure(1)
PyPlot.plot(t, z)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"z", rotation=0, fontsize=14)
PyPlot.title("Figure 1 Double pendulum with the second elastic: \$z\$ against t", fontsize=16)
PyPlot.savefig("graphics/Figure 1 Double pendulum with the second elastic z vs t.png")
if @isdefined(p2)
    PyPlot.close(p2)
end
p2 = PyPlot.figure(2)
PyPlot.plot(t, theta1)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\theta_1", rotation=0, fontsize=14)
PyPlot.title("Figure 2 Double pendulum with the second elastic: \$\\theta_1\$ against t", fontsize=16)
PyPlot.savefig("graphics/Figure 2 Double pendulum with the second elastic theta1 vs t.png")
if @isdefined(p3)
    PyPlot.close(p3)
end
p3 = PyPlot.figure(3)
PyPlot.plot(t, theta2)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\theta_2", rotation=0, fontsize=14)
PyPlot.title("Figure 3 Double pendulum with the second elastic: \$\\theta_2\$ against t", fontsize=16)
PyPlot.savefig("graphics/Figure 3 Double pendulum with the second elastic theta2 vs t.png")
if @isdefined(p4)
    PyPlot.close(p4)
end
p4 = PyPlot.figure(4)
PyPlot.plot(t, dz)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\dot{z}", rotation=0, fontsize=14)
PyPlot.title("Figure 4 Double pendulum with the second elastic: \$\\dot{z}\$ against t", fontsize=16)
PyPlot.savefig("graphics/Figure 4 Double pendulum with the second elastic zdot vs t.png")
if @isdefined(p5)
    PyPlot.close(p5)
end
p5 = PyPlot.figure(5)
PyPlot.plot(t, dtheta1)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}_1", rotation=0, fontsize=14)
PyPlot.title("Figure 5 Double pendulum with the second elastic: \$\\dot{\\theta}_1\$ against t", fontsize=16)
PyPlot.savefig("graphics/Figure 5 Double pendulum with the second elastic theta1dot vs t.png")
if @isdefined(p6)
    PyPlot.close(p6)
end
p6 = PyPlot.figure(6)
PyPlot.plot(t, dtheta2)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}_2", rotation=0, fontsize=14)
PyPlot.title("Figure 6 Double pendulum with the second elastic: \$\\dot{\\theta}_2\$ against t", fontsize=16)
PyPlot.savefig("graphics/Figure 6 Double pendulum with the second elastic theta2dot vs t.png")
if @isdefined(p7)
    PyPlot.close(p7)
end
p7 = PyPlot.figure(7)
PyPlot.plot(z, dz)
PyPlot.xlabel(L"z", fontsize=14)
PyPlot.ylabel(L"\dot{z}", rotation=0, fontsize=14)
PyPlot.title("Figure 7 Double pendulum with the second elastic: \$\\dot{z}\$ against z", fontsize=16)
PyPlot.savefig("graphics/Figure 7 Double pendulum with the second elastic zdot vs z.png")
if @isdefined(p8)
    PyPlot.close(p8)
end
p8 = PyPlot.figure(8)
PyPlot.plot(theta2, dtheta2)
PyPlot.xlabel(L"\theta_2", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}_2", rotation=0, fontsize=14)
PyPlot.title("Figure 8 Double pendulum with the second elastic: \$\\dot{\\theta}_2\$ against \$\\theta_1\$", fontsize=16)
PyPlot.savefig("graphics/Figure 8 Double pendulum with the second elastic theta1dot vs theta1.png")
if @isdefined(p9)
    PyPlot.close(p9)
end
p9 = PyPlot.figure(9)
PyPlot.plot(theta2, dtheta2)
PyPlot.xlabel(L"\theta_2", fontsize=14)
PyPlot.ylabel(L"\dot{\theta}_2", rotation=0, fontsize=14)
PyPlot.title("Figure 9 Double pendulum with the second elastic: \$\\dot{\\theta}_2\$ against \$\\theta_2\$", fontsize=16)
PyPlot.savefig("graphics/Figure 9 Double pendulum with the second elastic theta2dot vs theta2.png")
include("double_pendulum_anim.jl")
double_pendulum_anim(t, params.r, params.l .+ z, theta1, theta2, 10001, 0.1, "Figure 10 Double pendulum 2nd elastic animation k is 1 without dissipation.mp4")