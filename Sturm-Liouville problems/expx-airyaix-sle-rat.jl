# My first Julia program
# Import package manager
using Pkg;

# Install and import LinearAlgebra
Pkg.add("LinearAlgebra")
using LinearAlgebra;

# Install and import SpecialFunctions
Pkg.add("SpecialFunctions")
using SpecialFunctions;

# Install and import PyPlot
Pkg.add("PyPlot")
using PyPlot;

N         = 1000;
prec      = 0.2435;
Nfrag     = Int(N*prec);
# Maximum index in ysub and the corresponding solution vector to be displayed in plot (as displaying the entire vectors makes it impossible to see any real details).
maxindex  = Int(3*N/5);
# 0.0438 for L=70; 0.1293 for L=60; 0.0446 for L=69;
L         = 1e2;
k         = 1;
# Column vector of integers from 0 to N
n         = 0:1:N;
# Chebyshev extrema grid
t         = pi*(n/N.+1);
x         = cos.(t);
xsub      = x[2:N];
y         = L*((cot.(t./2)).^2);
ysub      = y[2:N];
T         = cos.(acos.(x)*n');
Tsub      = T[2:N,:];
sqx       = -((xsub.^2).-1);
Usub      = Diagonal(((sqx).^(-0.5)))*sin.(acos.(xsub)*n');
dTsub     = Usub * Diagonal(n);
dT        = [(-((-1).^n).*n.^2)'; dTsub; (n.^2)'];
dT        = Diagonal(2*L./(((y.+L).^2)))*dT;
D1        = dT/T;
# Second-order differentiation matrix for extrema grid
D2        = D1*D1;
# Second-order differentiation matrix for extrema grid without endpoints
E1        = D1[2:N,2:N];
E2        = D2[2:N,2:N];
H         = E2 - 2*E1 - Diagonal(ysub.-1);
EIG       = eigen(H);
Y         = EIG.vectors;
Lam       = EIG.values;
Y         = Y[:, sortperm(real(Lam), rev=true)];
Y         = [zeros(1,N-1); Y; zeros(1,N-1)];
DY        = D1*Y;
Lam       = sort(real(Lam), rev=true);
Lamerr    = abs.(airyai.(Lam[1:Nfrag]));
Lamerrrms = sqrt(Lamerr'*Lamerr/Nfrag)

# Accuracy is rather lousy, only getting <10 eigenvalues to 10 decimal places with N=1e4. I have tried a few L values and it seems like the inaccuracy is here to stay. 