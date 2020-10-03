using PyPlot;
include("RKF45.jl");

# Parameter object
struct paramObj
    sigma::Float64
    beta::Float64
    rho::Float64
end

# Function representing the RHS of the problem
function f(params::paramObj, t::Float64, vars::Array{Float64, 1})
    sigma = params.sigma;
    beta  = params.beta;
    rho   = params.rho;
    x     = vars[1];
    y     = vars[2];
    z     = vars[3];
    dx    = sigma*(y-x);
    dy    = x*(rho-z) - y;
    dz    = x*y - beta*z;
    return [dx, dy, dz];
end

# Problem parameters
sigma = 10.0;
beta = 8.0/3.0;
rho = 28.0;
params = paramObj(sigma, beta, rho);

# Initial conditions and domain of integration
t0 = 0.0;
tf = 60.0;
x0 = 1.0;
y0 = 1.0;
z0 = 1.0;
conds = [x0 y0 z0];

# Error tolerance and initial step size
epsilon = 1e-9;
dtInitial = 0.1;

# Solve problem
solution = RKF45(f, params, t0, tf, conds, epsilon, dtInitial);

# Extract solution values
vars = solution.vars;
x = vars[:,1];
y = vars[:,2];
z = vars[:,3];

# Plot solution
PyPlot.figure(1)
PyPlot.clf();
PyPlot.plot3D(x, y, z);