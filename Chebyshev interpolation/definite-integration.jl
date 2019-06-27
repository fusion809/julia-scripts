using PyPlot

# Clear all variables
clear all()
# Show numbers in a long; yet intuitive form
format long g
# Number of collocation points
N         = 1000
NN        = 100000
# Integration interval
xx0       = 0
xxf       = 50
# Vector
n         = 0:N
# Extrema grid()
x         = -cos(pi*n'/N)
# Transformation to [xx0, xxf]
xtrans    = (xxf-xx0)/2*(x+1)
# linearly spaced grid()
xlin      = range(-1,1,length=NN+1)'
Tlin      = cos(acos(xlin)*n)
xtranslin = (xlin+1)*(xxf-xx0)/2
# Chebyshev polys of the first kind
T         = cos(acos(x)*n)
# Now for arrays that do not include the endpoints
xsub      = x[2:N]
Tsub      = T[2:N,:]
Usub      = diag(1./sqrt(1-xsub.^2))*sin(acos(xsub)*n)
dTsub     = Usub*diag(n)
# Add the endpoints
dT        = [-(n.^2).*(-1).^(n) dTsub  [n].^2]'
# LHS without a [coefficient vector]
H         = 2/(xxf-xx0)*dT
# RHS; essentially what we're trying to integrate
F         = gsimp[xtrans]
# Initial condition; with definite integrals y[xx0] = 0
H[1,:]    = T[1,:]
F[1]      = 0
# coefficient vector
a         = H\F
dya       = T\F
# solution to ODE vector
y         = T*a
# definite integral
y         = y-y[1]
# y on the linear grid()
ylin      = Tlin*a
dylin     = Tlin*dya
# exact indefinite solution per Wolfram Alpha
#yexact    = -1/10*exp(-xtranslin).*(-2*sin(2*xtranslin)+cos(2*xtranslin)+5)
#yexact    = real(yexact)
# exact definite solution
#yexact    = yexact-yexact[1]
# error to Chebyshev approximation
#err       = abs(yexact-ylin)
#errrms    = sqrt(err'*err/(NN+1))
# y approximated by quadratic integration function
#yquad     = quad("g", 0, inf)
# y solved from ODE
lsode_options["absolute tolerance", 1e-18]
yode      = lsode["g", 0, xtranslin]
# error in approximation
#errquad   = abs(yexact[end]-yquad)
diffodech = abs(yode-ylin)
errodech  = sqrt(diffodech'*diffodech/(NN+1))
figure(1)
semilogy(xtranslin,diffodech)
figure(2)
plot(xtranslin,ylin,"r',xtranslin[2:end],dylin[2:end],'g")
