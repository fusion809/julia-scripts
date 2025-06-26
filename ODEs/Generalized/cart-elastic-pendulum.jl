include("RKF45.jl")
params       = (g=9.81, l=1, k=1, m1=1, m2=1, b1=0.0, b2=0.0, c1=0.00, c2=0.00)
theta0       = pi/2;
dtheta0      = 0.0;
z0           = 0.0;
dz0          = 0.0;
s0           = 10.0;
ds0          = 0.0;
t0           = 0.0;
tf           = 100.0;
conds        = @SVector [z0, dz0, s0, ds0, theta0, dtheta0];
dtInit       = 0.01;
epsilon      = 1e-7;

function cartEP(params::NamedTuple, t::Float64, coords::SVector{6,Float64})::SVector{6,Float64}
    g        = params.g;
    l        = params.l;
    k        = params.k;
    m1       = params.m1;
    m2       = params.m2;
    b1       = params.b1;
    b2       = params.b2;
    c1       = params.c1;
    c2       = params.c2;
    z        = coords[1];
    dz       = coords[2];
    s        = coords[3];
    ds       = coords[4];
    theta    = coords[5];
    dtheta   = coords[6];
    vbob     = sqrt(ds^2+dz^2+(l+z)^2*dtheta^2+2*ds*(dz*sin(theta)+dtheta*(l+z)*cos(theta)));
    Qtheta   = -(b2 + c2*vbob)*(l+z)*(ds*cos(theta)+(l+z)*dtheta);
    Qs       = -(b1 + c1*abs(ds))*ds - (b2 + c2*vbob)*(ds + dz*sin(theta) + (l+z)*dtheta*cos(theta));
    Qz       = -(b2 + c2*vbob)*(ds*sin(theta) + dz);
    d2s      = 1/m1 *(Qs - Qtheta*cos(theta)/(l+z) - (k*z+Qz)*sin(theta));
    d2z      = (l+z)*dtheta^2 + g*cos(theta) - sin(theta)/m1*(Qs-Qtheta*cos(theta)/(l+z)) + (m1+m2*(sin(theta))^2)/(m1*m2)*(k*z+Qz);
    d2theta  = 1/(l+z)^2 * (Qtheta*(m1+m2*(cos(theta))^2)/(m1*m2) - (l+z)*(2*dz*dtheta + cos(theta)/m1*(Qs-(k*z+Qz)*sin(theta))+g*sin(theta)));
    return [dz, d2z, ds, d2s, dtheta, d2theta]
end

t, coords = RKF45(cartEP, params, t0, tf, conds, epsilon, dtInit, "relative");
using PyPlot;
if (@isdefined(p1))
    PyPlot.close(p1)
end
p1=PyPlot.figure(1);
PyPlot.plot(t, coords[:,1])
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"z")
if (@isdefined(p2))
    PyPlot.close(p2)
end
p2=PyPlot.figure(2);
PyPlot.plot(t, coords[:,3])
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"s")
if (@isdefined(p3))
    PyPlot.close(p3)
end
p3=PyPlot.figure(3);
PyPlot.plot(t, coords[:,5])
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"\theta")