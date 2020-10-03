"""
    f(t, x, y, z)

The right-hand side of the coupled system of 3 ODEs defined by the equation:

``\\dfrac{d\\vec{r}}{dt} = f(t, x, y, z)``.
"""
function f(delta, t, S, I, R)
    beta  = 1/3;
    gamma = 1/4;
    N     = 90;
    dS    = -beta*I*(1-delta)*S/N;
    dI    = beta*I*(1-delta)*S/N-gamma*I*(1-delta);
    dR    = gamma*I*(1-delta);
    return [dS, dI, dR];
end

"""
    RKF45(f::Function, dtInitial::Float64, epsilon::Float64, t0::Float64, tf::Float64, x0::Float64, y0::Float64, z0::Float64)

A function that, using the [Runge-Kutta-Fehlberg method](https://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method) with fourth order approximations and fifth order error checking, approximates the solution to the ODE system:

``\\dfrac{d\\vec{r}}{dt} = f(t, x, y, z)``

with adaptive step size. 
"""
function RKF45(f, dtInitial, delta, epsilon, t0, tf, S0, I0, R0)
    # Initialize globals
    i = 1;
    t = Float64[t0];
    S = Float64[S0];
    I = Float64[I0];
    R = Float64[R0];
    dt = dtInitial;

    # Loop steps in the domain
    while t[i] < tf
        # dt should be the smallest out of tf-t[i] and the dt determined
        # the last iteration, as otherwise we won't finish at exactly tf
        dt = min(dt, tf-t[i]);
        
        # r approximators
        K1 = dt*f.(delta, t[i], S[i], I[i], R[i]);
        k1 = K1[1];
        l1 = K1[2];
        m1 = K1[3];

        K2 = dt*f.(delta, t[i]+dt/4, S[i]+k1/4, I[i]+l1/4, R[i]+m1/4);
        k2 = K2[1];
        l2 = K2[2];
        m2 = K2[3];

        K3 = dt*f.(delta, t[i]+3*dt/8, S[i]+3*k1/32+9*k2/32, I[i]+3*l1/32+9*l2/32, R[i]+3*m1/32+9*m2/32);
        k3 = K3[1];
        l3 = K3[2];
        m3 = K3[3];

        K4 = dt*f.(delta, t[i]+12*dt/13, S[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, I[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, R[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197);
        k4 = K4[1];
        l4 = K4[2];
        m4 = K4[3];

        K5 = dt*f.(delta, t[i]+dt, S[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, I[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, R[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104);
        k5 = K5[1];
        l5 = K5[2];
        m5 = K5[3];

        K6 = dt*f.(delta, t[i]+dt/2, S[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, I[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, R[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40);
        k6 = K6[1];
        l6 = K6[2];
        m6 = K6[3];

        # x1, y1 and z1 are 4th order approximations to x, y, z
        # x2, y2 and z2 are 5th order approximations
        x1 = S[i] + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5;
        y1 = I[i] + 25*l1/216 + 1408*l3/2565 + 2197*l4/4104 - l5/5;
        z1 = R[i] + 25*m1/216 + 1408*m3/2565 + 2197*m4/4104 - m5/5;
        x2 = S[i] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55;
        y2 = I[i] + 16*l1/135 + 6656*l3/12825 + 28561*l4/56430 - 9*l5/50 + 2*l6/55;
        z2 = R[i] + 16*m1/135 + 6656*m3/12825 + 28561*m4/56430 - 9*m5/50 + 2*m6/55;

        # A representation of the error in our x approximation
        Rx = abs(x1-x2)/dt;
        Ry = abs(y1-y2)/dt;
        Rz = abs(z1-z2)/dt;
        RRKF45 = max(Rx, Ry, Rz);
        
        # What our step sixe should be multiplied by in order to achieve 
        # an error tolerance of epsilon
        sx = 0.84*(epsilon/Rx)^(1/4);
        sy = 0.84*(epsilon/Ry)^(1/4);
        sz = 0.84*(epsilon/Rz)^(1/4);
        s = min(sx, sy, sz);

        # If R is less than or equal to epsilon, move onto the next step, 
        # otherwise repeat the iteration with our corrected h value
        if RRKF45 <= epsilon
            push!(t, t[i] + dt);
            push!(S, x1);
            push!(I, y1);
            push!(R, z1);
            i += 1;
            dt *= s;
        else
            dt *= s;
        end
    end
    
    return [t, S, I, R]
end

# Initial conditions, bounds of integration, error tolerance and initial dt
# value
t0 = 0;
tf = 2000;
S0 = 89;
I0 = 1;
R0 = 0;
delta = 0.9;
epsilon = 1e-8;
dtInitial = 0.1;

# Solve the problem and write it to t, x, y, and z arrays
t, S, I, R = RKF45(f, dtInitial, delta, epsilon, t0, tf, S0, I0, R0);

# Plot the solution
using PyPlot
PyPlot.figure(1)
PyPlot.clf()
PyPlot.plot3D(S,I,R)
#PyPlot.figure(2)
#PyPlot.clf()
##PyPlot.plot(x,dy)
#PyPlot.figure(3)
#PyPlot.clf()
#PyPlot.plot(y,dy)
