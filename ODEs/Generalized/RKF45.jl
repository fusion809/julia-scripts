# StaticArrays help optimize this code
using StaticArrays;

"""
    RKF45(f::Function, params, t0::Number, tf::Number, conds::Vector, epsilon, dtInitial)

f should be a function that returns the RHS of the ODE being solved expressed as a system of first-order ODEs
in an array of the form `[dx1/dt, dx2/dt, dx3/dt, dx4/dt, ..., dxn/dt]`. Its arguments should be: params 
(an object containing problem parameters), t (a Float64) and `vars::Vector` (a column vector of the form 
`[element1; element2; element3; ...; elementn]`).

`params` should be a named tuple containing parameter values. For the simple
pendulum problem with pendulum length 1 metre, for example, it can be written
as:

`params = (g = 9.81, l = 1.0)`.

`t0` is the value of t (the independent variable) at the beginning of the integration.

`tf` is the value of t at the end of the integration.

`conds` is an SVector containing initial conditions. For the simple pendulum
problem, for example, the following code may be used (where theta0 and
thetaDot0) have been defined elsewhere:

`conds = @SVector [theta0, thetaDot0]`.

`epsilon` is the error tolerance for the problem.

`dtInitial` is the initial guess for dt.
"""
function RKF45(f::Function, params::NamedTuple, t0::Float64, tf::Float64, conds::SVector, epsilon::Float64, dtInitial::Float64, tolType::String = "absolute", dtMin::Float64 = (tf-t0)/1e8)
    # Initialize relevant variables
    dt = dtInitial;
    t = Float64[t0];
    vars = [conds];
    i = 1;
    ti = t0;

    # Loop over t under the solution for tf has been found
    while ( ti < tf )
        varsi =  vars[i];
        dt = minimum((dt, tf-ti));

        # RKF45 approximators
        K1 = dt*f(params, ti, varsi);
        K2 = dt*f(params, ti + dt/4, varsi + K1/4);
        K3 = dt*f(params, ti + 3*dt/8, varsi + 3*K1/32 + 9*K2/32);
        K4 = dt*f(params, ti + 12*dt/13, varsi + 1932*K1/2197 - 7200*K2/2197 + 7296*K3/2197);
        K5 = dt*f(params, ti + dt, varsi + 439*K1/216 - 8*K2 + 3680*K3/513 - 845*K4/4104);
        K6 = dt*f(params, ti + dt/2, varsi - 8*K1/27 + 2*K2 - 3544*K3/2565 + 1859*K4/4104 - 11*K5/40);

        # 4/5th order approximations to next step value
        vars1 = varsi + 25*K1/216 + 1408*K3/2565 + 2197*K4/4104 - K5/5;
        vars2 = varsi + 16*K1/135 + 6656*K3/12825 + 28561*K4/56430 - 9*K5/50 + 2*K6/55;

        # Determine if error is small enough to move on to next step
        if (tolType in ["relative", "rel", "R", "r", "Rel", "Relative"])
            R = maximum(abs.(vars2 - vars1)/abs.(vars1))/dt;
        elseif (tolType in ["absolute", "abs", "A", "a", "Abs", "Absolute"])
            R = maximum(abs.(vars2 - vars1))/dt;
        else
            error("tolType is set to an invalid value ($tolType), so exiting...")
        end
        s = (epsilon/(2*R))^(0.25);
        if (R <= epsilon)
            Base.push!(t, ti+dt);
            StaticArrays.push!(vars, vars1);
            i += 1;
            ti = t[i];
        end
        dt *= s;
        if (dt < dtMin)
            @warn("dt has reached $dt at t=$ti which is less than dtMin=$dtMin")
            if (tolType in ["absolute", "abs", "A", "a", "Abs", "Absolute"])
                tolType = "relative";
                @warn("As you are using an absolute tolerance type, we will switch to relative tolerance to see if this fixes the problem...")
            elseif (tolType in ["relative", "rel", "R", "r", "Rel", "Relative"])
                @warn("Breaking out of loop as tolerance type is already set to relative.")
                break
            else
                error("tolType is set to an invalid value ($tolType), so exiting...")
            end
        end
    end

    # Transpose and enter into NamedTuple
    vars = transpose(reduce(hcat, vars));
    return t, vars;
end