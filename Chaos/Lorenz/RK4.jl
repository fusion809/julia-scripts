using Pkg;

# Import PyPlots
Pkg.add("PyPlot")
using PyPlot;

N     = 10000000;
sigma = 10.0
rho   = 28.0
beta  = 8.0/3.0;
t0    = 0;
tf    = 80;
dt    = (tf-t0)/N;
t     = t0:dt:tf;

function f(beta, rho, sigma, t, x, y, z)

	# The Lorenz equations
	dx_dt = sigma.*(y - x);
	dy_dt = -y + x.*(-z + rho);
	dz_dt = - beta*z + x*y;

	# Return the derivatives as a vector
	return dx_dt, dy_dt, dz_dt
end;

function RK4(beta, rho, sigma, dt, t, x, y, z)
    K1 = dt.*f(beta, rho, sigma, t, x, y, z);
    k1 = K1[1];
    l1 = K1[2];
    m1 = K1[3];
    K2 = dt.*f(beta, rho, sigma, t + dt/2, x + k1/2, y + l1/2, z + m1/2);
    k2 = K2[1];
    l2 = K2[2];
    m2 = K2[3];
    K3 = dt.*f(beta, rho, sigma, t + dt/2, x + k2/2, y + l2/2, z + m2/2);
    k3 = K3[1];
    l3 = K3[2];
    m3 = K3[3];
    K4 = dt.*f(beta, rho, sigma, t + dt/2, x + k3, y + l3, z + m3);
    k4 = K4[1];
    l4 = K4[2];
    m4 = K4[3];
    dx = 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    dy = 1/6 * (l1 + 2*l2 + 2*l3 + l4);
    dz = 1/6 * (m1 + 2*m2 + 2*m3 + m4);
    
    return dx, dy, dz
end

x     = zeros(N+1);
x[1]  = -2.0;
y     = zeros(N+1);
y[1]  = 3.0;
z     = zeros(N+1);
z[1]  = 0.5;

for i=1:N
    diff = RK4(beta, rho, sigma, dt, t[i], x[i], y[i], z[i]);
    x[i+1] = x[i] + diff[1];
    y[i+1] = y[i] + diff[2];
    z[i+1] = z[i] + diff[3];
end


PyPlot.figure(1)
PyPlot.plot3D(x, y, z);
PyPlot.xlabel("x")
PyPlot.ylabel("y")
PyPlot.zlabel("z")
