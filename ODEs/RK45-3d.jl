"""
    f(t, x, y, z)

The right-hand side of the coupled system of 3 ODEs defined by the equation:

``\\dfrac{d\\vec{r}}{dt} = f(t, x, y, z)``.
"""
function f(t, x, y, z)
    sigma = 10;
    beta  = 8/3;
    rho   = 28;
    dx    = sigma*(y-x);
    dy    = x*(rho-z) - y;
    dz    = x*y - beta*z;
    return [dx, dy, dz];
end

"""
    RK45(f::Function, dtInitial::Float64, epsilon::Float64, t0::Float64, tf::Float64, x0::Float64, y0::Float64, z0::Float64)

A function that, using the [Runge-Kutta-Fehlberg method](https://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method) with fourth order approximations and fifth order error checking, approximates the solution to the ODE system:

``\\dfrac{d\\vec{r}}{dt} = f(t, x, y, z)``

with adaptive step size. 
"""
function RK45(f, dtInitial, epsilon, t0, tf, x0, y0, z0)
    i = 1;
    t = Float64[t0];
    x = Float64[x0];
    y = Float64[y0];
    z = Float64[z0];
    dt = dtInitial;
    while t[i] < tf
        # dt should be the smallest out of tf-t[i] and the dt determined
        # the last iteration, as otherwise we won't finish at exactly tf
        dt = min(dt, tf-t[i]);
        
        # r approximators
        K1 = dt*f.(t[i], x[i], y[i], z[i]);
        k1 = K1[1];
        l1 = K1[2];
        m1 = K1[3];

        K2 = dt*f.(t[i]+dt/4, x[i]+k1/4, y[i]+l1/4, z[i]+m1/4);
        k2 = K2[1];
        l2 = K2[2];
        m2 = K2[3];

        K3 = dt*f.(t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32, z[i]+3*m1/32+9*m2/32);
        k3 = K3[1];
        l3 = K3[2];
        m3 = K3[3];

        K4 = dt*f.(t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, z[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197);
        k4 = K4[1];
        l4 = K4[2];
        m4 = K4[3];

        K5 = dt*f.(t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, z[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104);
        k5 = K5[1];
        l5 = K5[2];
        m5 = K5[3];

        K6 = dt*f.(t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, z[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40);
        k6 = K6[1];
        l6 = K6[2];
        m6 = K6[3];

        # x1, y1 and z1 are 4th order approximations to x, y, z
        # x2, y2 and z2 are 5th order approximations
        x1 = x[i] + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5;
        y1 = y[i] + 25*l1/216 + 1408*l3/2565 + 2197*l4/4104 - l5/5;
        z1 = z[i] + 25*m1/216 + 1408*m3/2565 + 2197*m4/4104 - m5/5;
        x2 = x[i] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55;
        y2 = y[i] + 16*l1/135 + 6656*l3/12825 + 28561*l4/56430 - 9*l5/50 + 2*l6/55;
        z2 = z[i] + 16*m1/135 + 6656*m3/12825 + 28561*m4/56430 - 9*m5/50 + 2*m6/55;

        # A representation of the error in our x approximation
        Rx = abs(x1-x2)/dt;
        Ry = abs(y1-y2)/dt;
        Rz = abs(z1-z2)/dt;
        R = max(Rx, Ry, Rz);
        
        # What our step sixe should be multiplied by in order to achieve 
        # an error tolerance of epsilon
        sx = 0.84*(epsilon/Rx)^(1/4);
        sy = 0.84*(epsilon/Ry)^(1/4);
        sz = 0.84*(epsilon/Rz)^(1/4);
        s = min(sx, sy, sz);

        # If R is less than or equal to epsilon, move onto the next step, 
        # otherwise repeat the iteration with our corrected h value
        if R <= epsilon
            push!(t, t[i] + dt);
            push!(x, x1);
            push!(y, y1);
            push!(z, z1);
            i += 1;
            dt *= s;
        else
            dt *= s;
        end
    end
    
    return [t, x, y, z]
end

# Initial conditions, bounds of integration, error tolerance and initial dt
# value
t0 = 0;
tf = 200;
x0 = 10;
y0 = 10;
z0 = 10;
epsilon = 1e-8;
dtInitial = 0.1;
t, x, y, z = RK45(f, dtInitial, epsilon, t0, tf, x0, y0, z0);
using PyPlot
PyPlot.figure(1)
PyPlot.clf()
PyPlot.plot3D(x,y,z)
#PyPlot.figure(2)
#PyPlot.clf()
##PyPlot.plot(x,dy)
#PyPlot.figure(3)
#PyPlot.clf()
#PyPlot.plot(y,dy)