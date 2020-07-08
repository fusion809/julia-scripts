function f(t,x, y)
    alpha = 1.1; beta=0.4;
    delta = 0.1; gamma = 0.4;
    return [alpha*x-beta*x*y, delta*x*y-gamma*y];
end

function rk45(t0, tf, x0, y0)
    epsilon = 1e-10;
    h = 0.1;
    t = Float64[t0];
    tfinal = tf;
    x = Float64[x0];
    y = Float64[y0];
    i = 1;
    while t[i]<tfinal
        h = min(h, tfinal-t[i]);
        K1 = h*f(t[i],x[i], y[i]);
        k1 = K1[1];
        l1 = K1[2];
        K2 = h*f(t[i]+h/4, x[i]+k1/4, y[i]+l1/4);
        k2 = K2[1];
        l2 = K2[2];
        K3 = h*f(t[i]+3*h/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32);
        k3 = K3[1];
        l3 = K3[2];
        K4 = h*f(t[i]+12*h/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197);
        k4 = K4[1];
        l4 = K4[2];
        K5 = h*f(t[i]+h, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104);
        k5 = K5[1];
        l5 = K5[2];
        K6 = h*f(t[i]+h/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40);
        k6 = K6[1];
        l6 = K6[2];
        x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        y1 = y[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        y2 = y[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
        R = abs(x1-x2)/h;
        s = 0.84*(epsilon/R)^(1/4);
        if R<=epsilon
            push!(t, t[i]+h)
            push!(x, x2)
            push!(y, y2)
            i = i+1;
            h = s*h;
        else
            h = s*h;
        end
    end
    return [t, x, y]
end

t0 = 0;
tf = 30;
x0 = 10;
y0 = 10;
t = nothing; x = nothing; y = nothing; xx = nothing; yy=nothing; splx=nothing; sply=nothing; tt = nothing;
t, x, y = rk45(t0, tf, x0, y0);
tt = LinRange(t0, tf, Int64(1e7+1));
using Pkg;
Pkg.add("Dierckx")
using Dierckx
splx = Spline1D(t, x);
xx = splx(tt);
sply = Spline1D(t, y);
yy = sply(tt);
using PyPlot
PyPlot.figure(1)
PyPlot.clf()
PyPlot.plot(tt,xx, :red, label="Number of prey animals")
PyPlot.plot(tt,yy, :blue, label="Number of predator animals")
PyPlot.xlabel("Time")
PyPlot.ylabel("Species population")
PyPlot.xlim((-0.1,30.1))
PyPlot.ylim((-0.1,31.5))
PyPlot.legend(loc="upper right",fancybox="true")
PyPlot.figure(2)
PyPlot.clf()
PyPlot.plot(xx,yy)
PyPlot.xlabel("Number of prey animals")
PyPlot.ylabel("Number of predator animals")