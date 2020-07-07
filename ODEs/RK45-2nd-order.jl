using FunctionIntegrator

function f(x,y, dy)
    return [dy, -9.8*cos(y)];
end

function rk45(t0, tf, theta0, dtheta0)
    epsilon = 1e-12;
    h = 0.1;
    x = Float64[t0];
    xfinal = tf;
    y = Float64[theta0];
    dy = Float64[dtheta0];
    error_in_y = Float64[0];
    i = 1;
    while x[i]<xfinal
        h = min(h, xfinal-x[i]);
        K1 = h*f(x[i],y[i], dy[i]);
        k1 = K1[1];
        l1 = K1[2];
        K2 = h*f(x[i]+h/4, y[i]+k1/4, dy[i]+l1/4);
        k2 = K2[1];
        l2 = K2[2];
        K3 = h*f(x[i]+3*h/8, y[i]+3*k1/32+9*k2/32, dy[i]+3*l1/32+9*l2/32);
        k3 = K3[1];
        l3 = K3[2];
        K4 = h*f(x[i]+12*h/13, y[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, dy[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197);
        k4 = K4[1];
        l4 = K4[2];
        K5 = h*f(x[i]+h, y[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, dy[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104);
        k5 = K5[1];
        l5 = K5[2];
        K6 = h*f(x[i]+h/2, y[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, dy[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40);
        k6 = K6[1];
        l6 = K6[2];
        y1 = y[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        dy1 = dy[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        y2 = y[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        dy2 = dy[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
        R = abs(y1-y2)/h;
        s = 0.84*(epsilon/R)^(1/4);
        if R<=epsilon
            push!(x, x[i]+h)
            push!(y, y1)
            push!(dy, dy1)
            i = i+1;
            h = s*h;
        else
            h = s*h;
        end
    end
    return [x, y, dy]
end

t0 = 0;
tf = 2*chebyshev_quadrature(x -> 1/sqrt(-19.6*sin(x)), 20, 1, -pi, 0);
theta0 = 0;
dtheta0 = 0;
x, y, dy = rk45(t0, tf, theta0, dtheta0);
xx = LinRange(t0, tf, Int64(1e7+1));
using Pkg;
Pkg.add("Dierckx")
using Dierckx
spl = Spline1D(x, y);
yy = spl(xx);
error = findmin(yy)[1]+pi;
using PyPlot
PyPlot.figure(1)
PyPlot.clf()
PyPlot.plot(xx,yy)
PyPlot.figure(2)
PyPlot.clf()
PyPlot.plot(x,dy)
PyPlot.figure(3)
PyPlot.clf()
PyPlot.plot(y,dy)