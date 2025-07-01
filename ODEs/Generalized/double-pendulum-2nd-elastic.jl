include("RKF45.jl")
params = (m1 = 1, m2=1, g=9.81, r=1, k=1, l=1, 
b1=0.0, c1=0.00, b2=0.0, c2=0.00);
z0 = 0.0;
dz0 = 0.0;
theta10 = 0;
dtheta10 = 0;
theta20 = 0;
dtheta20 = 0;
t0 = 0.0;
tf = 200.0;
conds = @SVector [theta10, dtheta10, theta20, dtheta20, z0, dz0];
dtInit = 1e-3;
epsilon = 1e-5;

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
    theta1   = coords[1];
    dtheta1  = coords[2];
    theta2   = coords[3];
    dtheta2  = coords[4];
    z        = coords[5];
    dz       = coords[6];
    v1       = r*abs(dtheta1);
    v2       = sqrt(r^2*dtheta1^2 + dz^2 + (l+z)^2*dtheta2^2 + 2*r*dtheta1 * dz * sin(theta2-theta1) + 2*r*(l+z)*dtheta1*dtheta2*cos(theta2-theta1));
    bob1dr   = -(b1 + c1*v1);
    bob2dr   = -(b2 + c2*v2);
    Qz       = bob2dr * (r*dtheta1 * sin(theta2-theta1) + dz);
    Qtheta1  = bob1dr * (r^2*dtheta1) + bob2dr * (r^2*dtheta1 + r*dz*sin(theta2-theta1) + r*(l+z)*dtheta2*cos(theta2-theta1));
    Qtheta2  = bob2dr * (l+z) * (r*dtheta1*cos(theta2-theta1) + (l+z)*dtheta2);
    firstout = -(m1+m2*(cos(theta2-theta1))^2)/m1;
    d2theta2 = firstout * (Qtheta2/(m2*(l+z)^2)-g*cos(theta2)/(l+z)) - (r*dtheta1^2*sin(theta2-theta1)+2*dz*dtheta2)/(l+z) - (cos(theta2-theta1))/(m1*r*(l+z)) * (m2*g*r*sin(theta2)*sin(theta2-theta1) + Qtheta1 + r*(k*z-Qz)*sin(theta2-theta1)-(m1+m2)*g*r*cos(theta1));
    outer    = 1/((m1+m2*(cos(theta2-theta1))^2)*r^2);
    inner    = Qtheta1 - g*(m1+m2)*r*cos(theta1) + r*(k*z-Qz)*sin(theta2-theta1) - m2*r*((r*dtheta1^2*cos(theta2-theta1))*sin(theta2-theta1) + (2*dz*dtheta2 + (l+z)*dtheta2)*cos(theta2-theta1));
    d2theta1 = outer * inner;
    d2z      = -r*d2theta1 * sin(theta2-theta1) + r*dtheta1^2*cos(theta1-theta1) + (l+z)*dtheta2^2 - g*sin(theta2) + (Qz-k*z)/m2;
    return [dtheta1, d2theta1, dtheta2, d2theta2, dz, d2z];
end

t, coords = RKF45(DP2E, params, t0, tf, conds, epsilon, dtInit);
theta1    = coords[:,1];
dtheta1   = coords[:,2];
theta2    = coords[:,3];
dtheta2   = coords[:,4];
z         = coords[:,5];
dz        = coords[:,6];
using PyPlot
if @isdefined(p1)
    PyPlot.close(p1)
end
p1 = PyPlot.figure(1)
PyPlot.plot(t, theta1)
PyPlot.xlabel(L"t", fontsize=14)
PyPlot.ylabel(L"\theta_1", rotation=0, fontsize=14)
PyPlot.title("Double pendulum with the second elastic: \$\\theta_1\$ against t", fontsize=16)