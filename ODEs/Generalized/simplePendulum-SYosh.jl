const t0 = 0.0;
const tf = 60.0;
const dt = 0.0001; 
using PyPlot;
const m = 1.0;
const r = 1.0;
const g = 9.81;
const theta0 = 0;
const p0 = 0; 
include("SYosh.jl")

function fp(params,theta)
    m = params.m;
    g = params.g;
    r = params.r;
    return -m*g*r*cos.(theta)
end

function fq(params,p)
    m = params.m;
    r = params.r;
    return p./(m*r^2)
end

params = (g=g, r=r, m=m);
conds = @SVector [theta0, p0];
t, theta, p = SYosh(fp, fq, params, t0, tf, conds, dt);

using PyPlot;
PyPlot.figure(1)
PyPlot.plot(t, theta)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"\theta", rotation=0)
PyPlot.title("Simple pendulum: integrated using Yoshida's 4th-order symplectic method")
PyPlot.figure(2)
PyPlot.plot(t, p)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"p", rotation=0)
PyPlot.title("Simple pendulum: integrated using Yoshida's 4th-order symplectic method")
PyPlot.figure(3)
PyPlot.plot(theta, p)
PyPlot.xlabel(L"\theta")
PyPlot.ylabel(L"p", rotation=0)
PyPlot.title("Simple pendulum: integrated using Yoshida's 4th-order symplectic method")
H = p.^2 ./ (2 * m * r^2) .+ m * g * r .* sin.(theta)
PyPlot.figure(4)
PyPlot.plot(t, H)
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"H", rotation=0)
PyPlot.title("Hamiltonian over time")