# My first Julia program
N         = 1000;
a         = 0;
b         = 100;
k         = 1;

n         = linspace(0,N,N+1)';
nn        = collect(0:N);
x         = -1 + 2*n'/N;
xsub      = x[2:N];
ysub      = (b-a)/2*xsub + (a+b)/2;
T         = cos(acos(x)*n);
Tsub      = T[2:N,:];
Usub      = *(diagm(1./sqrt(1-xsub.^2)), sin(acos(xsub)*n));
dTsub     = *(Usub, diagm(nn));
dT        = [-((-1).^n).*n.^2; dTsub; n.^2];
# Second derivative of T(x) on extrema grid without endpoints
d2Tsubbr  = *(diagm(xsub), Usub) - *(Tsub, diagm(nn));
d2Tsubb   = *(diagm(1./(1-xsub.^2)), d2Tsubbr);
d2Tsub    = *(d2Tsubb, diagm(nn));
# Second derivative of T(x) on extrema grid with endpoints
d2T       = [((-1).^n).*(n.^2).*(n.^2-1)/3; d2Tsub; (n.^2).*(n.^2-1)/3];
# Second-order differentiation matrix for extrema grid
D2        = *(d2T, inv(T));
# Second-order differentiation matrix for extrema grid without endpoints
E2        = D2[2:N,2:N];
H         = -4/((b-a)^2) * E2 + k * diagm(ysub);
LAM, Y    = eig(H);
LAM
