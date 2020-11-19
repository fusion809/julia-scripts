using PyPlot;
include("RKF45.jl");

# Function representing the RHS of the problem
function Lorenz(params::NamedTuple, t::Float64, 
    vars::SVector{3,Float64})::SVector{3,Float64}
    sigma, beta, rho = params.sigma, params.beta, params.rho;
    x, y, z = vars[1], vars[2], vars[3];
    dx    = sigma*(y-x);
    dy    = x*(rho-z) - y;
    dz    = x*y - beta*z;
    return [dx, dy, dz];
end

# Problem parameters
params = (sigma = 10.0, beta = 8.0/3.0, rho = 28.0);

# Initial conditions and domain of integration
t0 = 0.0;
tf = 200.0;
x0, y0, z0 = 1.0, 1.0, 1.0;
conds = @SVector [x0, y0, z0];

# Error tolerance and initial step size
epsilon = 1e-9;
dtInitial = 0.1;

# Solve problem
@time begin
t, vars = RKF45(Lorenz, params, t0, tf, conds, epsilon, dtInitial);
end

# Extract solution values
x, y, z = vars[:,1], vars[:,2], vars[:,3];

# Plot solution
PyPlot.figure(1);
PyPlot.clf();
PyPlot.plot3D(x, y, z);
PyPlot.title("Lorenz attractor");
PyPlot.figure(2);
PyPlot.clf();
PyPlot.plot(t, x, label=L"x");
PyPlot.plot(t, y, label=L"y");
PyPlot.plot(t, z, label=L"z");
PyPlot.legend();