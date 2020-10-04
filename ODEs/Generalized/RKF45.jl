# Solution object
struct solObj
    t::Array{Float64,1}
    vars::Array{Float64,2}
end

"""
    RKF45(f::Function, params, t0::Number, tf::Number, conds::Vector, epsilon::Float64, dtInitial::Float64)

f should be a function that returns the RHS of the ODE being solved expressed as a system of first-order ODEs
in an array of the form [dx1/dt, dx2/dt, dx3/dt, dx4/dt, ..., dxn/dt]. Its arguments should be: params 
(an object containing problem parameters), t (a Float64) and vars::Vector (a column vector of the form 
[element1; element2; element3; ...; elementn]).

params should be an object containing parameter values to be fed to f. It should be created using:

```julia
struct paramObj
    param1
    param2
    param3
    ...
    paramn
end
params = paramObj(param1, param2, param3, ..., paramn);
```

t0 is the value of t (the independent variable) at the beginning of the integration.
tf is the value of t at the end of the integration.
conds is a row vector of the form [x1(0) x2(0) x3(0) x4(0) ... xn(0)].
epsilon is the error tolerance for the problem.
dtInitial is the initial choice for dt.
"""
function RKF45(f::Function, params::NamedTuple, t0::Float64, tf::Float64, conds::Array{Float64}, epsilon::Float64, dtInitial::Float64)
    # Initialize global variables
    dt = dtInitial;
    t = Float64[t0];
    vars = [conds];
    i = 1;
    ti = t0;

    # Loop through each time step
    while ( ti < tf )
        varsi =  vars[i]
        dt = minimum((dt, tf-ti))
        K1 = dt*f(params, ti, varsi)
        K2 = dt*f(params, ti + dt/4, varsi + K1/4)
        K3 = dt*f(params, ti + 3*dt/8, varsi + 3*K1/32 + 9*K2/32)
        K4 = dt*f(params, ti + 12*dt/13, varsi + 1932*K1/2197 - 7200*K2/2197 + 7296*K3/2197)
        K5 = dt*f(params, ti + dt, varsi + 439*K1/216 - 8*K2 + 3680*K3/513 - 845*K4/4104)
        K6 = dt*f(params, ti + dt/2, varsi - 8*K1/27 + 2*K2 - 3544*K3/2565 + 1859*K4/4104 - 11*K5/40)
        vars1 = varsi + 25*K1/216 + 1408*K3/2565 + 2197*K4/4104 - K5/5
        vars2 = varsi + 16*K1/135 + 6656*K3/12825 + 28561*K4/56430 - 9*K5/50 + 2*K6/55
        R = maximum(abs.(vars2 - vars1))/dt
        s = (epsilon/(2*R))^(0.25)
        if (R <= epsilon)
            push!(t, ti+dt)
            push!(vars, vars1)
            i += 1
            ti = t[i]
        end
        dt *= s
    end
    vars = transpose(reduce(hcat, vars));
    return solObj(t, vars);
end