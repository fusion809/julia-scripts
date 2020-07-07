# Elastic pendulum
function f(t, x, dx, theta, dtheta)
    g  = 9.8;
    l0 = 1.0;
    m  = 1.0;
    k  = 1.0;
    d2x = (l0+x)*dtheta^2-k*x/m+g*cos(theta);
    d2theta = -(g/(l0+x))*sin(theta) - 2*dx*dtheta/(l0+x);
    return [dx, d2x, dtheta, d2theta];
end

function rk45(h, i, t, x, dx, theta, dtheta)
    K1 = h*f.(t[i], x[i], dx[i], theta[i], dtheta[i]);
    k1 = K1[1];
    l1 = K1[2];
    m1 = K1[3];
    n1 = K1[4];

    K2 = h*f.(t[i]+h/4, x[i]+k1/4, dx[i]+l1/4, theta[i]+m1/4, dtheta[i]+n1/4);
    k2 = K2[1];
    l2 = K2[2];
    m2 = K2[3];
    n2 = K2[4];

    K3 = h*f.(t[i]+3*h/8, x[i]+3*k1/32+9*k2/32, dx[i]+3*l1/32+9*l2/32, theta[i]+3*m1/32+9*m2/32, dtheta[i]+3*n1/32+9*n2/32);
    k3 = K3[1];
    l3 = K3[2];
    m3 = K3[3];
    n3 = K3[4];

    K4 = h*f.(t[i]+12*h/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, dx[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, theta[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197, dtheta[i]+1932*n1/2197-7200*n2/2197+7296*n3/2197);
    k4 = K4[1];
    l4 = K4[2];
    m4 = K4[3];
    n4 = K4[4];

    K5 = h*f.(t[i]+h, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, dx[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, theta[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104, dtheta[i]+439*n1/216-8*n2+3680*n3/513-845*n4/4104);
    k5 = K5[1];
    l5 = K5[2];
    m5 = K5[3];
    n5 = K5[4];

    K6 = h*f.(t[i]+h/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, dx[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, theta[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40, dtheta[i]-8*n1/27+2*n2-3544*n3/2565+1859*n4/4104-11*n5/40);
    k6 = K6[1];
    l6 = K6[2];
    m6 = K6[3];
    n6 = K6[4];

    x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
    dx1 = dx[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
    theta1  = theta[i] + 25*m1/216+1408*m3/2565+2197*m4/4104-m5/5;
    dtheta1 = dtheta[i] + 25*n1/216+1408*n3/2565+2197*n4/4104-n5/5;
    x2      = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
    dx2     = dx[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
    theta2  = theta[i] + 16*m1/135+6656*m3/12825+28561*m4/56430-9*m5/50+2*m6/55;
    dtheta2 = dtheta[i] + 16*n1/135+6656*n3/12825+28561*n4/56430-9*n5/50+2*n6/55;

    R = abs(x1-x2)/h;
    return R, x2, dx2, theta2, dtheta2
end

function solver(t0, tf, x0, dx0, theta0, dtheta0)
    epsilon = 1e-9;
    h       = 0.1;
    t       = Float64[t0];
    x       = Float64[x0];
    dx      = Float64[dx0];
    theta   = Float64[theta0];
    dtheta  = Float64[dtheta0]
    i       = 1;
    while t[i] < tf
        h = min(h, tf-t[i]);
        R, x2, dx2, theta2, dtheta2 = rk45(h, i, t, x, dx, theta, dtheta);
        s = 0.84*(epsilon/R)^(1/4);

        if R<=epsilon
            push!(t, t[i]+h)
            push!(x, x2)
            push!(dx, dx2)
            push!(theta, theta2)
            push!(dtheta, dtheta2)
            i = i+1;
        end
        h = s*h;
    end
    return [t, x, dx, theta, dtheta]
end

t0 = 0;
tf = 100;
x0 = 1;
dx0 = 0;
theta0 = 3*pi/4;
dtheta0 = 0;
t, x, dx, theta, dtheta = solver(t0, tf, x0, dx0, theta0, dtheta0);
using PyPlot
PyPlot.figure(1)
PyPlot.clf()
PyPlot.plot3D(t, x, theta)
#PyPlot.figure(2)
#PyPlot.clf()
##PyPlot.plot(x,dy)
#PyPlot.figure(3)
#PyPlot.clf()
#PyPlot.plot(y,dy)