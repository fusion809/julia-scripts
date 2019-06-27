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

N         = 10000;
# Maximum index in ysub and the corresponding solution vector to be displayed in plot (as displaying the entire vectors makes it impossible to see any real details).
maxindex  = Int(3*N/5);
prec      = 0.27;
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
# @ L=366, prec=0.27, eigerrrms=6.70e-11
# @ L=367, prec=0.27, eigerrrms=1.05e-10
L         = 366;
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
H         = - E2 + k * Diagonal(ysub);
EIG       = eigen(H);
Y         = EIG.vectors;
Lam       = EIG.values;
Y         = Y[:, sortperm(Lam)];
Lam       = sort(Lam);
eigerr    = abs.(airyai.(-Lam[1:Nfrag]));
eigerrrms = sqrt(eigerr'*eigerr/(Nfrag));
Y         = [zeros(1,N-1); Y; zeros(1,N-1)];
D1IY2     = D1\(Y.^2);
IntY2     = abs.(D1IY2[N+1,:]-D1IY2[1,:]);
normcoef  = (IntY2).^(-0.5);
Y         = Y*Diagonal(normcoef);

PyPlot.figure(1)
PyPlot.plot(ysub[1:maxindex],Y[1:maxindex,20])
PyPlot.figure(2)
PyPlot.semilogy(eigerr)
