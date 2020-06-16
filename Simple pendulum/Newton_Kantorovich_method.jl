# This script uses the Newton-Kantorovich method to approximate
# the solution to the problem of the simple pendulum.
# theta double dot = -g/l cos(theta)
# theta(0) = theta dot (0) = 0
# Our first guess is theta = pi/2*(cos(2*pi*t/T)+1)
# Then we iteratively solve the linear ODE:
# delta_i double dot - g/l sin(theta_i) delta_i = -(theta_i double dot +
# g/l*cos(theta_i))
# delta_i(0) = delta_i dot(0) = 0
# theta_(i+1) = theta_i + delta_i
# for i=1,2,3,...,NN
using Pkg;

# Import PyPlots
Pkg.add("PyPlot")
using PyPlot;

# Required to determine the period of the curve
Pkg.add("QuadGK")
using QuadGK;

# Installing required Julia modules
Pkg.add("LinearAlgebra")
using LinearAlgebra;

# N+1 is the number of points on our extrema grid
N                            = 100;
# NN iterations are used when applying the Newton-Kantorovich method
NN                           = 4;
# Acceleration due to gravity in metres per second squared
g                            = 9.8;
# length of the pendulum in metres
l                            = 1.0;
function f(y)
    # abs is used below to prevent complex number errors due to the sqrt
    return 1.0 / sqrt(abs(2.0 * g / l * sin(y)))
end
# The period of the simple pendulum, that is, how long it takes for
# it to go from theta=0 back to theta=0
period                       = 2*quadgk(y -> f(y), -pi, 0, rtol=1e-14)[1];
# Column vector of integers from 0 to N
n                            = 0:1:N;
# Chebyshev extrema grid; tt is the parameterization
tt                           = pi*(-n/N.+1.0);
# The grid itself
x                            = cos.(tt);
# Our transformed domain variable
t                            = period/2.0*(x.+1);
# Clearing x as it's unused henceforth
x                            = nothing;
# T_{mn} = T_n(x_m), where T_n is the Chebyshev polynomial of
# the first kind, and x_m are points on the Chebyshev extrema grid
T                            = cos.(tt*n');
# T'_n(x_m) = nU_{n-1}(x_m), m=1,2,3,...,N-1
dTsub                        = Diagonal(csc.(tt[2:N]))*sin.(tt[2:N]*n')*Diagonal(n);
# Adding endpoints, m=0 & N.
dT                           = [(-((-1.0).^n).*n.^2)'; dTsub; (n.^2)'];
# Clear unused variables to free up RAM
dTsub                        = nothing;
tt                           = nothing;
# Calculate first order differentiation matrix
D1                           = dT/T;
D1                           = 2.0/period*D1;
D2                           = D1^2;
theta                        = zeros(N+1,NN+1);
# A more precise approximation that we provide as we know what the curve looks
# like
theta[:,1]                   = pi/2*(cos.(2.0*pi*t/period).-1);
# A more improvised approximation, obtained using:
# theta ~ At^2
# cos(theta) ~ 1-theta^2/2
# And setting the residual at t=period/2 to zero
#theta[:,1]                   = -32/(period^4)*(l/g+sqrt(l^2/g^2+(period^4)/16))*t.^2;
# Negative residue is required as part of the Newton-Kantorovich method
negative_residue             = zeros(N+1,NN+1);

# Solving our linear ODEs for delta_i
for i=1:1:NN
    # Our Newton-Kantorovich linear ODE differential operator
    # Our linear ODE can be written as L delta = f
    # Where L is our differential operator and f is the negative residue
    differential_operator          = (D2-g/l*Diagonal(sin.(theta[:,i])));
    negative_residue[:,i]          = (-(D2*theta[:,i]+g/l*(cos.(theta[:,i]))));
    # Initial conditions; delta must meet them two in order for theta to
    differential_operator[1,1]     = 1;
    differential_operator[N+1,:]   = D1[1,:];
    negative_residue[1,i]          = 0;
    negative_residue[N+1,i]        = 0;
    # Adding our delta to theta[:,i] to get our new theta (theta[:,i+1]) values
    theta[:,i+1]                   = theta[:,i]+differential_operator\negative_residue[:,i];
end

# NNN+1 is how many points are in our linear grid
NNN                          = 100000;
# Basis function coefficients
a                            = T\theta[:,NN+1];
# Our linear grid we're using to get less pixelated plots
linear_grid                  = LinRange(-1,1,NNN+1);
transformed_linear_grid      = period/2.0*(linear_grid.+1);
# T_n(t) on our linear grid
T_linear_grid                = cos.(acos.(linear_grid)*n');
# theta on a linear grid
theta_linear_grid            = T_linear_grid*a;

# How much our iteration has improved our initial estimate of theta
correction                   = abs.(theta[:,NN+1]-theta[:,1]);
rms_correction               = sqrt(correction'*correction/(N+1));
# Substitute our values of theta on the extrema grid into the
# original equation and find the residue
residual                     = D2*theta+g/l*(cos.(theta));
rms_residual                 = zeros(NN+1,1);
# Root mean square of these values
for i=1:NN+1
    rms_residual[i]          = sqrt(residual[:,i]'*residual[:,i]/(N+1));
end

# Plots
# Enable separate graph windows
pygui(true)
# Compare our first guess of theta with our final solution
fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(t,theta[:,1])
ax2.plot(transformed_linear_grid,theta_linear_grid)
#figure(2)
#clf()
#semilogy(t,residual[:,1])
#figure(3)
#clf()
#semilogy(t,residual[:,6])
#figure(4)
#clf()
#semilogy(t,residual[:,NN+1])
figure(5)
clf()
semilogy(rms_residual)
xlabel(L"$i$ in $\theta_i$")
ylabel("RMS (residual)")
title("Semilog plot of the root mean square of the residual")
