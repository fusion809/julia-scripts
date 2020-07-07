function f(t, x, y, z)
    sigma = 10;
    beta  = 8/3;
    rho   = 28;
    dx    = sigma*(y-x);
    dy    = x*(rho-z) - y;
    dz    = x*y - beta*z;
    return [dx, dy, dz];
end

function rk45(h, i, t, x, y, z)
    K1 = h*f.(t[i], x[i], y[i], z[i]);
    k1 = K1[1];
    l1 = K1[2];
    m1 = K1[3];

    K2 = h*f.(t[i]+h/4, x[i]+k1/4, y[i]+l1/4, z[i]+m1/4);
    k2 = K2[1];
    l2 = K2[2];
    m2 = K2[3];

    K3 = h*f.(t[i]+3*h/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32, z[i]+3*m1/32+9*m2/32);
    k3 = K3[1];
    l3 = K3[2];
    m3 = K3[3];

    K4 = h*f.(t[i]+12*h/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, z[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197);
    k4 = K4[1];
    l4 = K4[2];
    m4 = K4[3];

    K5 = h*f.(t[i]+h, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, z[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104);
    k5 = K5[1];
    l5 = K5[2];
    m5 = K5[3];

    K6 = h*f.(t[i]+h/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, z[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40);
    k6 = K6[1];
    l6 = K6[2];
    m6 = K6[3];

    x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
    y1 = y[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
    z1 = z[i] + 25*m1/216+1408*m3/2565+2197*m4/4104-m5/5;
    x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
    y2 = y[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
    z2 = z[i] + 16*m1/135+6656*m3/12825+28561*m4/56430-9*m5/50+2*m6/55;
    
    R = abs(x1-x2)/h;
    return R, x2, y2, z2
end

function solver(t0, tf, x0, y0, z0)
    epsilon = 1e-9;
    h = 0.1;
    t = Float64[t0];
    x = Float64[x0];
    y = Float64[y0];
    z = Float64[z0];
    i = 1;
    while t[i] < tf
        h = min(h, tf-t[i]);
        R, x2, y2, z2 = rk45(h, i, t, x, y, z);
        s = 0.84*(epsilon/R)^(1/4);

        if R<=epsilon
            push!(t, t[i]+h)
            push!(x, x2)
            push!(y, y2)
            push!(z, z2)
            i = i+1;
        end
        h = s*h;
    end
    return [t, x, y, z]
end

t0 = 0;
tf = 100;
x0 = 1;
y0 = 1;
z0 = 1;
t, x, y, z = solver(t0, tf, x0, y0, z0);
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