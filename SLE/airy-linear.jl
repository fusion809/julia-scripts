# Approximate the solution to the problem
# -d2y/dx2 + kxy = lambda y
# using Chebyshev spectral methods with a linear
# transformation from the Chebyshev extrema grid

# Import package manager
using Pkg;

# Installing required Julia modules
Pkg.add("LinearAlgebra")
using LinearAlgebra;

# Required to use the Airy function and its derivative
Pkg.add("SpecialFunctions")
using SpecialFunctions;

# Matplotlib plotting function
Pkg.add("PyPlot");
Pkg.add("PyCall");
using PyPlot;
pygui(true)

# Shouldn't need installing, but needed for useful plot annotations
using Printf;

# Our Chebyshev spectral grid will have N+1 points
N                       = 10000;
# The number of eigenvalues and eigenvectors we're computing
Nfrag                   = 5361;
# Our truncated integration domain is [a,b]
a                       = 0.0;
b                       = 870.0;
# The coefficient in -d2y/dx2 + kxy = lambda y
k                       = 1.0;
# Column vector of integers from 0 to N
n                       = 0:1:N;
# Chebyshev extrema grid
t                       = pi*(-n/N.+1.0);
x                       = cos.(t);
# Our transformed domain variable
y                       = (b-a)/2*x.+(a+b)/2;
# Clearing x as it's unused henceforth
x                       = nothing;
# T_{mn} = T_n(x_m), where T_n is the Chebyshev polynomial of
# the first kind, and x_m are points on the Chebyshev extrema grid
T                       = cos.(t*n');
# T'_n(x_m) = nU_{n-1}(x_m), m=1,2,3,...,N-1
dTsub                   = Diagonal(csc.(t[2:N]))*sin.(t[2:N]*n')*Diagonal(n);
# Adding endpoints, m=0 & N.
dT                      = [(-((-1.0).^n).*n.^2)'; dTsub; (n.^2)'];
dTsub                   = nothing;
# Calculate first order differentiation matrix
D1                      = dT/T;
# Second order differentiation matrix
D2                      = D1*D1;
D1                      = nothing;
# Second-order differentiation matrix for extrema grid without endpoints
E2                      = D2[2:N,2:N];
# Clearing D2
D2                      = nothing;
# Hy = lambda y is our problem
H                       = -4/((b-a)^2) * E2 + k * Diagonal(y[2:N]);
# Clearing E2
E2                      = nothing;
# Eigenfunctions and eigenvalues for the operator H
EIG                     = eigen(H);
# Clear H to free up RAM
H                       = nothing;
# Extract eigenfunctions and eigenvalues from EIG
eigenfunction_approx    = EIG.vectors;
eigenvalue_approx       = EIG.values;
# Clear EIG to free up RAM
EIG                     = nothing;
# Sort both in order of ascending absolute value of eigenvalues
IX                      = sortperm(eigenvalue_approx, by=abs);
eigenfunction_approx    = eigenfunction_approx[:, IX];
eigenvalue_approx       = eigenvalue_approx[IX];
IX                      = nothing;
# Add endpoints to eigenfunction approximations
eigenfunction_approx    = [zeros(1,N-1); eigenfunction_approx; zeros(1,N-1)];

# Use Newton's method to refine our eigenvalues
function newtons(xinput)
    Ai      = airyai(-xinput);
    Aip     = airyaiprime(-xinput);
    xoutput = xinput;
    while abs(Ai/Aip) > 1e-12
        Ai      = airyai(-xoutput);
        Aip     = airyaiprime(-xoutput);
        xoutput = xoutput + Ai/Aip;
    end
    Ai  = nothing;
    Aip = nothing;
    return xoutput
end

# eigenvalues that have been refined through Newton's method
eigenvalue_exact           = eigenvalue_approx[1:Nfrag];
# eigenfunctions calculated using the analytical solution and
# eigenvalue_exact
eigenfunction_exact          = zeros(N+1,Nfrag);
# root mean square of the error in eigenfunctions
rms_of_eigenfunction_error   = zeros(N+1,1);
# error in eigenfunctions
eigenfunction_error          = zeros(N+1,Nfrag);
eigenfunction_scaling_factor = zeros(N+1,1);
# loop is required to calculate some things
for i in 1:1:Nfrag
    # Refine eigenvalue approximations using Newton's method
    eigenvalue_exact[i]             = newtons(eigenvalue_exact[i]);
    # An eigenfunction approximation that uses the analytical solution
    eigenfunction_exact[:,i]        = airyai.(y.-eigenvalue_exact[i])
    # Our Chebyshev-approximated eigenfunctions may be off from our
    # analytical ones by a constant multiplier, so let's scale it
    eigenfunction_scaling_factor[i] = (sign(eigenfunction_exact[2,i])/sign(eigenfunction_approx[2,i]))*(findmax(abs.(eigenfunction_exact[:,i]))[1]/findmax(abs.(eigenfunction_approx[:,i]))[1])
    eigenfunction_approx[:,i]       = eigenfunction_scaling_factor[i]*eigenfunction_approx[:,i];
    # Relative error in eigenfunctions (absolute error/
    # max of exact eigenfunction)
    eigenfunction_error[:,i]        = abs.(eigenfunction_exact[:,i]-eigenfunction_approx[:,i])./findmax(eigenfunction_exact[:,i])[1];
    # Root mean square of eigenfunction relative error
    rms_of_eigenfunction_error[i,1] = sqrt(eigenfunction_error[:,i]'*eigenfunction_error[:,i]/(N+1));
end
# Clear eigenfunction_error to save RAM
eigenfunction_error          = nothing;
# Clear eigenfunction_scaling_factor to save RAM
eigenfunction_scaling_factor = nothing;
# Error in the Chebyshev approximation of the eigenvalues
eigenvalue_error             = abs.((eigenvalue_approx[1:Nfrag]-eigenvalue_exact)./eigenvalue_exact);
# Calculate root mean square of eigenvalue relative error
rms_of_eigenvalue_error      = sqrt(eigenvalue_error'*eigenvalue_error/(Nfrag));
# The root mean square of the root mean square of the error in
# eigenfunctions
rms2_of_eigenfunction_error  = sqrt(rms_of_eigenfunction_error'*rms_of_eigenfunction_error/(N+1))[1];

# Print this root mean square error in eigenvalues
print("RMS of eigenvalue error is ", rms_of_eigenvalue_error, "\n")
print("RMS of RMS of eigenfunction error is ", rms2_of_eigenfunction_error, "\n")

# Plots; setting it up as a function as it makes it easier to
# call it again in the notebook if required.
function eigenplots()
    # Eigenfunction plots
    PyPlot.figure(1)
    PyPlot.clf()
    PyPlot.xlim((y[1],y[801]))
    PyPlot.plot(y[1:801],eigenfunction_approx[1:801,1])
    PyPlot.title(latexstring("Plot of the ", L"1^\mathrm{st}", " eigenfunction, corresponding to ", L"\lambda = ",
    @sprintf("%.10g", eigenvalue_approx[1])))
    PyPlot.figure(2)
    PyPlot.clf()
    PyPlot.xlim((y[1],y[1201]))
    PyPlot.plot(y[1:1201],eigenfunction_approx[1:1201,10])
    PyPlot.title(latexstring("Plot of the ", L"10^\mathrm{th}", " eigenfunction, corresponding to ", L"\lambda = ",
    @sprintf("%.10g", eigenvalue_approx[10])))
    PyPlot.figure(3)
    PyPlot.clf()
    PyPlot.xlim((y[1],y[2501]))
    PyPlot.plot(y[1:2501],eigenfunction_approx[1:2501,100])
    PyPlot.title(latexstring("Plot of the ", L"100^\mathrm{th}", " eigenfunction, corresponding to ", L"\lambda = ",
    @sprintf("%.10g", eigenvalue_approx[100])))
    PyPlot.figure(4)
    PyPlot.clf()
    PyPlot.xlim((y[1],y[2501]))
    PyPlot.plot(y[1:2501],eigenfunction_approx[1:2501,200])
    PyPlot.title(latexstring("Plot of the ", L"200^\mathrm{th}", " eigenfunction, corresponding to ", L"\lambda = ",
    @sprintf("%.10g", eigenvalue_approx[200])))
    # Semilog plot of eigenvalue errors
    PyPlot.figure(5)
    PyPlot.clf()
    PyPlot.xlim((1,Nfrag))
    PyPlot.semilogy(eigenvalue_error)
    PyPlot.title("Semilog plot of eigenvalue errors")
    # Semilog plot of root mean square of eigenfunction errors
    PyPlot.figure(6)
    PyPlot.clf()
    # The following line is required, as otherwise
    # x from 0 to 10,000 is shown.
    PyPlot.xlim((1,Nfrag))
    PyPlot.semilogy(rms_of_eigenfunction_error)
    PyPlot.title("Semilog plot of the root mean square of eigenfunction error")
end

# Commented out, as plots are no longer needed
eigenplots()
