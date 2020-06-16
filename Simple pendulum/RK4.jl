# This script uses the 4th order Runge-Kutta method to
# numerically approximate the solution to the ODE:
# \ddot{\theta} = - \dfrac{g}{l} \cos{\theta}
# \theta(0) = \dot{\theta}(0) = 0
# Install and import necessary packages
using Pkg;
Pkg.add("QuadGK")
using QuadGK;
Pkg.add("PyPlot")
using PyPlot;
pygui(true)

# Number of steps
N = 10000000;
# Acceleration due to gravity in metres per second squared
g = 9.8;
# Length of the pendulum in metres
l = 1.0;

# t = \int_{\theta_0}^{\theta_1} h(\theta) d\theta
function h(x)
    return -1/sqrt(abs(2*g/l*sin(x)))
end

# Determine the period of the problem
T      = 2*quadgk(x -> h(x), 0, -pi, rtol=1e-14)[1];

# Integration domain is [t0,tf]
t0     = 0;
tf     = T;
dt     = (tf-t0)/N;
t      = t0:dt:tf;

# Initiate our solution vectors
theta  = zeros(N+1,1);
# dtheta = \dot{\theta}
dtheta = zeros(N+1,1);

# Loop to apply RK4
for i=1:N
    # k1, k2, k3, k4 are used to approximate theta
    # l1, l2, l3, l4 are for dtheta
    k1 = dt*dtheta[i];
    l1 = dt*(-g/l*cos(theta[i]));
    k2 = dt*(dtheta[i]+l1/2);
    l2 = dt*(-g/l*cos(theta[i]+k1/2));
    k3 = dt*(dtheta[i]+l2/2);
    l3 = dt*(-g/l*cos(theta[i]+k2/2));
    k4 = dt*(dtheta[i]+l3);
    l4 = dt*(-g/l*cos(theta[i]+k3));
    # Determine next theta and dtheta values
    theta[i+1] = theta[i] + 1/6*(k1+2*k2+2*k3+k4);
    dtheta[i+1] = dtheta[i] + 1/6*(l1+2*l2+2*l3+l4);
end

error_theta_min = abs(findmin(theta)[1]+pi);
error_dtheta_min = abs(findmin(dtheta)[1] + sqrt(2*g/l));

# We know from integrating our original equation with
# respect to theta that:
# \dot{\theta} = \pm \sqrt{-2g/l \sin{\theta}}
# this simply calculates the difference between our
# RK4 calculated theta dot and this analytical
# expression for it based on our RK4-calculated
# theta values.
residual_dtheta = abs.(abs.(dtheta)-abs.(sqrt.(abs.(2*g/l*sin.(theta)))));
rms_residual_dtheta = sqrt(residual_dtheta'*residual_dtheta/(N+1))[1];

# Plots
figure(1)
clf()
plot(t,theta)
figure(2)
clf()
plot(t,dtheta)
# Phase plot
figure(3)
clf()
plot(theta,dtheta)
# semilog plot of residual in dtheta
figure(4)
clf()
semilogy(t,residual_dtheta)

println("error_theta_min is $(error_theta_min).")
println("error_dtheta_min is $(error_dtheta_min).")
println("rms_residual_dtheta is $(rms_residual_dtheta).")
