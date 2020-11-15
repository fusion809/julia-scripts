using PyPlot;
include("RKF45.jl");

# Function representing the RHS of the ODE system
function SIR(params::NamedTuple, t::Float64, vars::SVector{3,Float64})::SVector{3,Float64}
    S     = vars[1];
    I     = vars[2];
    R     = vars[3];
    beta  = params.beta;
    gamma = params.gamma;
    delta = params.delta;
    dS    = -beta*I*(1-delta)*S/N;
    dI    = beta*I*(1-delta)*S/N-gamma*I;
    dR    = gamma*I;
    return [dS, dI, dR];
end

# Define problem parameters
beta = 1.5;
gamma = 0.25;
delta = 0.9;
N = 100;
params = nothing;
params = (beta = beta, gamma = gamma, delta = delta, N = N);

# Define initial conditions and domain of integration
t0 = 0.0;
tf = 1e2;
S0 = 89.0;
I0 = 11.0;
R0 = 0.0;
conds = @SVector [S0, I0, R0];

# Step size and error tolerance
epsilon = 1e-12;
dtInitial = 0.1;

# Solve problem
@time begin
t, vars = RKF45(SIR, params, t0, tf, conds, epsilon, dtInitial);
end

# Extract solution values
S = vars[:,1];
I = vars[:,2];
R = vars[:,3];

# Plot solution
# Phase plot is not used, as for SIR it's pretty boring
PyPlot.figure(1)
PyPlot.clf();
PyPlot.plot(t, S);
PyPlot.plot(t, I);
PyPlot.plot(t, R);
