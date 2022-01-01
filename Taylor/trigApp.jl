# Use Taylor series to approximate trig functions
function cosApp(N, x)
    n = 0:N
    x = BigFloat.(x)
    ser = BigFloat.((-1).^n.*x.^(2*n))./BigFloat.(factorial.(big.(2*n)))
    app = sum(ser)
    err = app - cos(x)
    return app, err
end

function sinApp(N, x)
    n = 0:N
    x = BigFloat.(x);
    ser = BigFloat.((-1).^n.*x.^(2*n.+1))./BigFloat.(factorial.(big.(2*n.+1)))
    app = sum(ser)
    err = app - sin.(x)
    return app, err
end

using PyPlot
err = ones((30, 1))
for i in 1:length(err)
    app, err[i] = sinApp(i, pi)
end
PyPlot.figure(1)
PyPlot.semilogy(err)
PyPlot.show()