# My first Julia program
N         = 1000;
n         = linspace(0,N,N+1)';
x         = -1 + 2*n'/N;
xsub      = x[2:N]
T         = cos(acos(x)*n);
Tsub      = T[2:N,:]
Usub      = diagm(1./sqrt(1-xsub.^2))*sin(acos(xsub)*n);
dTsub     = Usub*diagm(n);
dT        = [-((-1).^n).*n.^2; dTsub; n.^2];
# Second derivative of T(x) on extrema grid without endpoints
d2Tsub    = diagm(1./(1-xsub.^2))*(diagm(xsub)*Usub-Tsub*diagm(n))*diagm(n);
# Second derivative of T(x) on extrema grid with endpoints
d2T       = [((-1).^n).*(n.^2).*(n.^2-1)/3; d2Tsub; (n.^2).*(n.^2-1)/3];
# Second-order differentiation matrix for extrema grid
D2        = d2T/T;
# Second-order differentiation matrix for extrema grid without endpoints
E2        = D2[2:N,2:N];
