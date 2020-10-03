# Solution object
struct solObj
    t::Vector
    vars::Array{Float64,2}
end

"""
    RKF45(f::Function, params, t0::Number, tf::Number, conds::Vector, epsilon::Float64, dtInitial::Float64)

f should be a function that returns the RHS of the ODE being solved expressed as a system of first-order ODEs
in an array of the form [dx1/dt, dx2/dt, dx3/dt, dx4/dt, ..., dxn/dt]. Its arguments should be: params 
(an object containing problem parameters), t (a Float64) and vars::Vector (a row vector of the form 
[element1 element2 element3 ... elementn]).

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
function RKF45(f::Function, params, t0::Float64, tf::Float64, conds::Array{Float64, 2}, epsilon::Float64, dtInitial::Float64)
    dt = dtInitial;
    t = Float64[t0];
    vars = conds;
    i = 1;
    while ( t[i] < tf )
        varsi = vars[i,:];
        dt = min(dt, tf-t[i]);

        # Interpolants
        K1 = dt*f(params, t[i], varsi);
        K2 = dt*f(params, t[i] + dt/4.0, varsi + K1/4.0);
        K3 = dt*f(params, t[i] + 3.0*dt/8.0, varsi + 3.0*K1/32.0 + 9.0*K2/32.0);
        K4 = dt*f(params, t[i] + 12*dt/13, varsi + 1932*K1/2197 - 7200*K2/2197 + 7296*K3/2197);
        K5 = dt*f(params, t[i] + dt, varsi + 439*K1/216 - 8*K2 + 3680*K3/513 - 845*K4/4104);
        K6 = dt*f(params, t[i] + dt/2, varsi - 8*K1/27 + 2*K2 - 3544*K3/2565 + 1859*K4/4104 - 11*K5/40);

        # 4th/5th order approximations
        vars1 = varsi + 25*K1/216 + 1408*K3/2565 + 2197*K4/4104 - K5/5;
        vars2 = varsi + 16*K1/135 + 6656*K3/12825 + 28561*K4/56430 - 9*K5/50 + 2*K6/55;

        # Adjust step size if needed
        R = maximum(abs.(vars2-vars1))/dt;
        s = (epsilon/(2*R))^(0.25);
        if ( R <= epsilon )
            t = vcat(t, t[i]+dt);
            vars = vcat(vars, vars1');
            i += 1;
        end
        dt *= s;
    end
    solution = solObj(t, vars);
    return solution;
end