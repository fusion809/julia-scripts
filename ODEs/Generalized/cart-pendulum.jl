include("RKF45.jl")
using DifferentialEquations
params = (m1 = 1, m2r=1, m2b=1, l=1, g=9.81, alpha=5*pi/6, bcart = 0.01, 
ccart = 0.001, brod = 0.01, crod = 0.001, bbob = 0.01, cbob = 0.001);
z0 = 100.0;
dz0 = 100.0;
theta0 = -pi;
dtheta0 = 2;
t0 = 0.0;
tf = 20.0;
conds = @SVector [z0, dz0, theta0, dtheta0];
dtInit = 1e-3;
epsilon = 1e-5;

function cartPen(params::NamedTuple, t::Float64, coords::SVector{4,Float64})::SVector{4,Float64}
    m1      = params.m1;
    m2r     = params.m2r;
    m2b     = params.m2b;
    l       = params.l;
    g       = params.g;
    alpha   = params.alpha;
    bcart   = params.bcart;
    ccart   = params.ccart;
    brod    = params.brod;
    crod    = params.crod;
    bbob    = params.bbob;
    cbob    = params.cbob; 
    z       = coords[1];
    dz      = coords[2];
    theta   = coords[3];
    dtheta  = coords[4];
    masssum = m1+m2r+m2b;
    massps1 = m2r + 2*m2b;
    massps2 = m2r + 3*m2b;
    # Friction-free
#    d2theta = massps1*cos(theta)/(2*(l^2*massps2/3+(massps1^2)*l*(sin(theta)^2)/(4*masssum)))*(g*cos(alpha)-massps1*l*(dtheta^2)*sin(theta)/(2*masssum));
#    d2z     = -g*sin(alpha) + massps1/(2*masssum)*l*(d2theta*sin(theta)+(dtheta^2)*cos(theta))
    denom = 4*m1*massps2 + (m2r^2)*(1+3*(cos(theta))^2+12*m2b^2*(1+(cos(theta))^2)+4*m2r*m2b*(2+3*(cos(theta))^2));
    vrod = sqrt(dz^2+(l^2*(dtheta)^2)/4+l*dz*dtheta*sin(theta));
    vbob = sqrt(dz^2+l^2*dtheta^2+2*l*dz*dtheta*sin(theta));
    # Rod coefficients
    rodco1 = dz + l*dtheta*sin(theta)/2;
    rodco2 = dz * sin(theta) + l*dtheta/2;
    # Bob coefficients
    bobco1 = dz + l*dtheta*sin(theta);
    bobco2 = dz * sin(theta) + l*dtheta;
    # Rod drags
    roddr1 = (brod + crod*vrod) * rodco1;
    roddr2 = (brod + crod*vrod) * rodco2;
    # Cart drag
    cartdr = (bcart + ccart*dz)*dz;
    # Bob drag
    bobdr1 = (bbob + cbob*vbob)*bobco1;
    bobdr2 = (bbob + cbob*vbob)*bobco2;
    firstterm = (6*g*massps1*masssum*cos(alpha)*cos(theta))/l;
    coef1 = 6*sin(theta)*massps1/l;
    secterm = massps1/2 * dtheta^2 * l*cos(theta);
    d2theta = 1/denom * (firstterm + coef1*(secterm-cartdr - roddr1 - bobdr1) - 6*masssum/l*(roddr2 + 2*bobdr2));
    d2z = -g*sin(alpha) - 1/masssum * (massps1*l/2 *(d2theta*sin(theta)+(dtheta^2)*cos(theta))-cartdr - roddr1 - bobdr1);
    return [dz, d2z, dtheta, d2theta];
end

function cartPen_wrapped!(u, p, t)
    return cartPen(p, t, u);
end

tspan = (t0, tf);
prob = ODEProblem{false}(cartPen_wrapped!, conds, tspan, params);
sol = solve(prob, Tsit5(), abstol=1e-11, reltol=1e-11);
t = sol.t;
u = sol.u;
z = zeros(length(u), 1);
dz = zeros(length(u), 1);
theta = zeros(length(u), 1);
dtheta = zeros(length(u), 1);
z      = getindex.(u, 1);
dz     = getindex.(u, 2);
theta  = getindex.(u, 3);
dtheta = getindex.(u, 4);
# t, coords = RKF45(cartPen, params, t0, tf, conds, epsilon, dtInit);
# using PyPlot;
# PyPlot.figure(1)
# PyPlot.plot(t, coords[:,1])
# PyPlot.xlabel(L"t")
# PyPlot.ylabel(L"z")
# PyPlot.figure(2)
# PyPlot.plot(t, coords[:,2])
# PyPlot.xlabel(L"t")
# PyPlot.ylabel(L"\dot{z}")
# PyPlot.figure(3)
# PyPlot.plot(t, coords[:,3])
# PyPlot.xlabel(L"t")
# PyPlot.ylabel(L"\theta")
# PyPlot.figure(4)
# PyPlot.plot(t, coords[:,4])
# PyPlot.xlabel(L"t")
# PyPlot.ylabel(L"\dot{\theta}")
using PyPlot;
function clearFigs(arr)
    for i in arr
        PyPlot.close(i)
    end
end
clearFigs([p1; p2; p3; p4])
p1 = PyPlot.figure(1);
PyPlot.plot(t, z)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"z")
p2 = PyPlot.figure(2);
PyPlot.plot(t, dz)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"\dot{z}")
p3 = PyPlot.figure(3);
PyPlot.plot(t, theta)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"\theta")
p4 = PyPlot.figure(4);
PyPlot.plot(t, dtheta)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"\dot{\theta}")



