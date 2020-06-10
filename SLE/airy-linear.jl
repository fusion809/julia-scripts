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
using PyPlot;
pygui(true)

# Shouldn't need installing, but needed for useful plot annotations
using Printf;

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
t                       = pi*(-n/N.+1);
x                       = cos.(t);
ysub                    = (b-a)/2*x[2:N].+(a+b)/2;
y                       = [a; ysub; b];
T                       = cos.(t*n');
# T'_n(x_m)
dT                      = [(-((-1).^n).*n.^2)'; Diagonal(((-((x[2:N].^2).-1)).^(-0.5)))*sin.(t[2:N]*n')*Diagonal(n); (n.^2)'];
# Calculate first order differentiation matrix
D1                      = dT/T;
# Clearing x as it's unused henceforth
x                       = nothing;
# Second order differentiation matrix
D2                      = D1*D1;
D1                      = nothing;
# Second-order differentiation matrix for extrema grid without endpoints
E2                      = D2[2:N,2:N];
# Clearing D2
D2                      = nothing;
# Hy = lambda y is our problem
H                       = -4/((b-a)^2) * E2 + k * Diagonal(ysub);
# Clearing ysub; unused in the rest of the script
ysub                    = nothing;
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
eigenfunction_approx    = eigenfunction_approx[:, sortperm(eigenvalue_approx, by=abs)];
eigenvalue_approx       = sort(eigenvalue_approx, by=abs);
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
eigenfunction_exact        = zeros(N+1,Nfrag);
# root mean square of the error in eigenfunctions
rms_of_eigenfunction_error = zeros(N+1,1);
# error in eigenfunctions
eigenfunction_error        = zeros(N+1,Nfrag);

# loop is required to calculate some things
for i in 1:1:Nfrag
    # Refine eigenvalue approximations using Newton's method
    eigenvalue_exact[i]             = newtons(eigenvalue_exact[i]);
    # An eigenfunction approximation that uses the analytical solution
    eigenfunction_exact[:,i]        = airyai.(y.-eigenvalue_exact[i])/abs(airyaiprime(-eigenvalue_exact[i]));
    # Our Chebyshev-approximated eigenfunctions may be off from our
    # analytical ones by a constant multiplier.
    eigenfunction_approx[:,i]       = eigenfunction_exact[2,i]/eigenfunction_approx[2,i]*eigenfunction_approx[:,i];
    # Relative error in eigenfunctions (absolute error/
    # max of exact eigenfunction)
    eigenfunction_error[:,i]        = abs.(eigenfunction_exact[:,i]-eigenfunction_approx[:,i])./findmax(eigenfunction_exact[:,i])[1];
    # Root mean square of eigenfunction relative error
    rms_of_eigenfunction_error[i,1] = sqrt(eigenfunction_error[:,i]'*eigenfunction_error[:,i]/(N+1));
end
# Clear eigenfunction_error to save RAM
eigenfunction_error     = nothing;
# Error in the Chebyshev approximation of the eigenvalues
eigenvalue_error        = abs.((eigenvalue_approx[1:Nfrag]-eigenvalue_exact)./eigenvalue_exact);
# Calculate root mean square of eigenvalue relative error
rms_of_eigenvalue_error = sqrt(eigenvalue_error'*eigenvalue_error/(Nfrag));

# Print this root mean square error in eigenvalues
print("RMS of eigenvalue error is ", rms_of_eigenvalue_error, "\n")

# Plots
PyPlot.figure(1)
PyPlot.xlim((y[1],y[801]))
PyPlot.plot(y[1:801],eigenfunction_approx[1:801,1])
PyPlot.title(latexstring("Plot of the ", L"1^\mathrm{st}", " eigenfunction, corresponding to ", L"\lambda = ", 
@sprintf("%.10f", eigenvalue_approx[1])))
PyPlot.figure(2)
PyPlot.xlim((y[1],y[1201]))
PyPlot.plot(y[1:1201],eigenfunction_approx[1:1201,10])
PyPlot.title(latexstring("Plot of the ", L"10^\mathrm{th}", " eigenfunction, corresponding to ", L"\lambda = ", 
@sprintf("%.10f", eigenvalue_approx[10])))
PyPlot.figure(3)
PyPlot.xlim((y[1],y[2501]))
PyPlot.plot(y[1:2501],eigenfunction_approx[1:2501,100])
PyPlot.title(latexstring("Plot of the ", L"100^\mathrm{th}", " eigenfunction, corresponding to ", L"\lambda = ", 
@sprintf("%.10f", eigenvalue_approx[100])))
PyPlot.figure(4)
PyPlot.xlim((y[1],y[2501]))
PyPlot.plot(y[1:2501],eigenfunction_approx[1:2501,200])
PyPlot.title(latexstring("Plot of the ", L"200^\mathrm{th}", " eigenfunction, corresponding to ", L"\lambda = ", 
@sprintf("%.10f", eigenvalue_approx[200])))
# Semilog plot of eigenvalue errors
PyPlot.figure(5)
PyPlot.xlim((1,Nfrag))
PyPlot.semilogy(eigenvalue_error)
PyPlot.title("Semilog plot of eigenvalue errors")
# Semilog plot of root mean square of eigenfunction errors
PyPlot.figure(6)
# The following line is required, as otherwise
# x from 0 to 10,000 is shown.
PyPlot.xlim((1,Nfrag))
PyPlot.semilogy(rms_of_eigenfunction_error)
PyPlot.title("Semilog plot of the root mean square of eigenfunction error")