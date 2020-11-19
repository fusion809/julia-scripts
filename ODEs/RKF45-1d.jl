function f(x,y)
    return 1+y^2
end

function rkf45()
    epsilon = 1e-13;
    h = 0.1;
    x = Float64[0.0];
    xfinal = pi/4;
    y = Float64[0];
    error_in_y = Float64[0];
    i = 1;
    while x[i]<xfinal
        h = min(h, xfinal-x[i]);
        k1 = h*f(x[i],y[i]);
        k2 = h*f(x[i]+h/4, y[i]+k1/4);
        k3 = h*f(x[i]+3*h/8, y[i]+3*k1/32+9*k2/32);
        k4 = h*f(x[i]+12*h/13, y[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197);
        k5 = h*f(x[i]+h, y[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104);
        k6 = h*f(x[i]+h/2, y[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40);
        y1 = y[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        y2 = y[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        R = abs(y1-y2)/h;
        s = 0.84*(epsilon/R)^(0.25);
        if R<=epsilon
            push!(x, x[i]+h)
            push!(y, y1)
            i += 1;
            push!(error_in_y, abs(y[i]-tan(x[i])));
        end
        h *= s;
    end
    return [x, y, error_in_y]
end

@time x, y, error = rkf45();
using PyPlot;
PyPlot.figure(1)
PyPlot.plot(x,y)
PyPlot.figure(2)
PyPlot.semilogy(error)
