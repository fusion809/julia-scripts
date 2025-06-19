const t0 = 0.0;
const tf = 300.0;
const dt = 0.0001; 
using PyPlot;
const m = 1.0;
const r = 1.0;
const g = 9.81;
const theta0 = 0.0;
const p0 = 0.0;
const k = 1:100; 
const period = 4*sqrt(r/g)*(pi/2 + 2*pi*sum(exp.(k*pi)./(1 .+ exp.(2*k*pi))))
include("SYosh.jl")
include("RKF45.jl")

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

function SP(params::NamedTuple, t::Float64, vars::SVector{2,Float64})::SVector{2,Float64}
    g = params.g;
    r = params.r;
    x = vars[1];
    dx = vars[2];
    array = [dx, -g/r*cos(x)];
    return array;
end

function Hamiltonian(params::NamedTuple, theta::Array{Float64}, p::Array{Float64})
    g = params.g;
    r = params.r;
    return p.^2 ./ (2 * m * r^2) .+ m * g * r .* sin.(theta);
end

function Hamiltonian_dth(params::NamedTuple, theta::Array{Float64}, dtheta::Array{Float64})
    g = params.g;
    r = params.r;
    p = m*r^2*dtheta;
    return Hamiltonian(params, theta, p);
end

params = (g=g, r=r, m=m);
conds = @SVector [theta0, p0];
epsilon = 1.5e-12;
t1, theta1, p1 = SYosh(fp, fq, params, t0, tf, conds, dt);
t2, TH = RKF45(SP, params, t0, tf, conds, epsilon, dt);
theta2 = TH[:, 1];
dtheta2 = TH[:, 2];

using PyPlot;
test = t1 .<= period;
PyPlot.figure(1)
PyPlot.plot(t1[test], theta1[test])
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"\theta", rotation=0)
PyPlot.title("Simple pendulum: integrated using Yoshida (1T)")
PyPlot.figure(2)
PyPlot.plot(t1[test], p1[test])
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"p", rotation=0)
PyPlot.title("Simple pendulum: integrated using Yoshida (1T)")
PyPlot.figure(3)
PyPlot.plot(theta1, p1)
PyPlot.xlabel(L"\theta")
PyPlot.ylabel(L"p", rotation=0)
PyPlot.title("Simple pendulum: integrated using Yoshida")
H1 = Hamiltonian(params, theta1, p1);
H2 = Hamiltonian_dth(params, theta2, dtheta2);
PyPlot.figure(4)
PyPlot.semilogy(t1, abs.(H1), label="Yoshida (symplectic)", color=:red)
PyPlot.semilogy(t2, abs.(H2), label="RKF45 (non-symplectic)", color=:blue)
PyPlot.legend()
PyPlot.xlabel(L"t")
PyPlot.ylabel(L"H", rotation=0)
PyPlot.title("Absolute value of the Hamiltonian over time")