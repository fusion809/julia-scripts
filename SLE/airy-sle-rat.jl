# Import package manager
using Pkg;

# Install and import LinearAlgebra
Pkg.add("LinearAlgebra")
using LinearAlgebra;

# Install and import SpecialFunctions
Pkg.add("SpecialFunctions")
using SpecialFunctions;

# The docs suggested that the following would cause Julia plots
# To be displayed in a Qt GUI, but it doesn't seem to work for 
# Jupyter Lab, or at least not if Jupyter Lab is run in a chroot.
#Pkg.add("PyCall")
#using PyCall;
#pygui(:qt)
# Install and import PyPlot
Pkg.add("PyPlot")
using PyPlot;

# Number of steps
N         = 10000;
# Maximum index in ysub and the corresponding solution vector to 
# be displayed in plot (as displaying the entire vectors makes 
# it impossible to see any real details).
maxindex  = Int(3*N/5);
# What percentage of the total eigenvalues and eigenvectors we
# wish to examine
prec      = 0.27;
# The corresponding number of eigenvalues and eigenvectors
Nfrag     = Int(N*prec);
# N=10,000
# At L=300; eigerrrms=4.56e-8 at prec=0.27
# At L=320; eigerrrms=1.73e-9 at prec=0.27
# At L=330; eigerrrms=4.82e-10 at prec=0.27
# At L=340; eigerrrms=2.45e-10 at prec=0.27
# At L=350; eigerrrms=1.12e-10 at prec=0.27
# At L=360, prec=0.27 eigerrrms=8.49e-11
# At L=370, prec=0.27 eigerrrms=1.06e-10
# @ L=365, prec=0.27, eigerrrms=7.89e-11
# @ L=365.9, prec=0.27, eigerrms=8.517802623658321e-11
# @ L=365.99, prec=0.27, eigerrrms=6.483350564293549e-11
# @ L=366, prec=0.27, eigerrrms=6.70e-11
# @ L=366.1, prec=0.27, eigerrrms=1.0010701888393075e-10
# @ L=367, prec=0.27, eigerrrms=1.05e-10
# With L=360, but our new relative Lamerrrms:
# prec=0.27, Lamerrrms = 1.726331851521917e-13
L         = 366;
k         = 1;
# Column vector of integers from 0 to N
n         = 0:1:N;
# Chebyshev extrema grid; 
# Do not use pi*(n/N.+1) as it causes
# our dT definition to fail
t         = pi*(-n/N.+1);
x         = cos.(t);
# Semi-infinite domain rational function transformation
y         = L*((cot.(t./2)).^2);
T         = cos.(t*n');
# do not switch acos(xsub) in the sin function below to t, 
# it gives far higher errors like 5 orders of magnitude
dT        = [(-((-1).^n).*n.^2)'; Diagonal(((-((x[2:N].^2).-1)).^(-0.5)))*sin.(t[2:N]*n')*Diagonal(n); (n.^2)'];
dT        = Diagonal(2*L./(((y.+L).^2)))*dT;
# Second-order differentiation matrix for extrema grid
D1        = dT/T;
D2        = D1*D1;
dT        = nothing;
T         = nothing;
# Second-order differentiation matrix for extrema grid without endpoints
#E1        = D1[2:N,2:N];
#E2        = D2[2:N,2:N];
EIG       = eigen(-D2[2:N,2:N] + k * Diagonal(y[2:N]));
D2        = nothing;
# Eigenvectors which are discretized versions of the 
# solution eigenfunctions
Y         = EIG.vectors;
# Eigenvalues
Lam       = EIG.values;
EIG       = nothing;
# Sort both by ascending absolute value of eigenvalues
Y         = Y[:, sortperm(Lam, by=abs)];
Lam       = sort(Lam, by=abs);

# The following lines of code that have been commented out
# pertain to a less precise way of estimating how inaccurate the
# eigenvalues are, and a precise way of estimating the error in the
# eigenfunctions.
#eigerr    = abs.(airyai.(-Lam[1:Nfrag]));
#eigerrrms = sqrt(eigerr'*eigerr/(Nfrag));
#Y         = [zeros(1,N-1); Y; zeros(1,N-1)];
#D1IY2     = D1\(Y.^2); # Integral of Y
#Y         = Y*Diagonal((abs.(D1IY2[N+1,:]-D1IY2[1,:])).^(-0.5));
#Yexact    = zeros(maxindex,Int(N/5));
#errrms    = zeros(maxindex,1);
#err       = zeros(maxindex,Int(N/5));
#Ymaxindex = Y[1:maxindex,:];
#for i in 1:1:1296
#    Yexact[:,i] = airyai.(y[1:maxindex].-Lam[i])/abs(airyaiprime(-Lam[i]));
#    err[:,i]=abs.(abs.(Yexact[:,i])-abs.(Ymaxindex[:,i]));
#    errrms[i,1]=sqrt(err[:,i]'*err[:,i]/maxindex);
#end

# Calculate more precise values of the eigenvalues
# We know from the analytical solution that the exact eigenvalues
# are the zeros of the function Ai(-x).
Lamex = Lam[1:Nfrag];
for i in 1:1:Nfrag
    Ai=1;
    while abs(Ai) > 1e-12
        Ai = airy(-Lamex[i]);
        Aip = airyprime(-Lamex[i]);
        Lamex[i] = Lamex[i] + Ai/Aip;
    end
end
Ai     = nothing;
Aip    = nothing;
Lamerr = abs.((Lam[1:Nfrag]-Lamex)./Lamex);
Lamerrrms = sqrt(Lamerr'*Lamerr/Nfrag);

Y         = [zeros(1,N-1); Y; zeros(1,N-1)];
D1IY2     = D1\(Y.^2); # Integral of Y
D1        = nothing;
Y         = Y*Diagonal((abs.(D1IY2[N+1,:]-D1IY2[1,:])).^(-0.5));
D1IY2     = nothing;
Yexact    = zeros(maxindex,Nfrag);
errrms    = zeros(Nfrag,1);
err       = zeros(maxindex,Nfrag);
for i in 1:1:Nfrag
    Yexact[:,i] = airyai.(y[1:maxindex].-Lamex[i])/abs(airyaiprime(-Lamex[i]));
    err[:,i]=abs.(abs.(Yexact[:,i])-abs.(Y[1:maxindex,i]));
    errrms[i,1]=sqrt(err[:,i]'*err[:,i]/maxindex);
end

print("Lamerrrms is ", Lamerrrms)
PyPlot.figure(1)
PyPlot.plot(y[1:maxindex],Y[1:maxindex,20])
PyPlot.figure(2)
PyPlot.semilogy(Lamerr)
plt.figure(3)
PyPlot.semilogy(errrms)
#PyPlot.figure(3)
#PyPlot.plot(y[1:maxindex],Yexact[:,20])
