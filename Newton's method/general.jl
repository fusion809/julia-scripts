#!/usr/bin/env julia
"""
    bisection(f::Function, N::Integer, a::Number, b::Number)

Uses the bisection method with N subdivisions of the specified domain [a, b] 
to find the roots of the function f within said domain. Used to find the 
initial guess that Newton's method then uses to get a more precise estimate of 
the roots.
"""
function bisection(f, N, a, b)
    x = LinRange(a, b, N+1)
    fval = f.(x)
    xv = x[2:end]
    initGuess = xv[diff(sign.(fval)).!= 0]
    return initGuess
end

"""
    newtons(f::Function, h::Float, tol::Float, itMax::Integer, 
    initGuess::Vector{Float})

Uses Newton's method (and if needed, up to third order Householder's methods) 
to refine initGuess, our initial guess of the root(s) of f, until either itMax 
iterations has been performed or the relative magnitude of the update Newton's 
provides to the estimate is below tol. h is the step size used to approximate 
the derivative.  
"""
function newtons(f, h, tol, itMax, initGuess)
    # Derivatives
    function fd(x, h)
        return (f(x+h)-f(x-h))/(2*h)
    end
    function f2d(x, h)
        return (f(x+h)-2*f(x)+f(x-h))/(h^2)
    end
    function f3d(x, h)
        return (f(x+2*h)-2*f(x+h)+f(x-h)-f(x-2*h))/(2*h^3)
    end

    # Check to see if (relative) difference in update is below tol
    function tolCheck(diff, sol, tol)
        if ! (sol ≈ 0)
            check = abs(diff/sol) > tol
        else
            check = abs(diff) > tol
        end
        return check
    end

    # Initialize sol
    sol = initGuess
    count = Int64.(zeros(size(initGuess)))
    for j=1:length(initGuess)
        diff = 1
        while (tolCheck(diff, sol[j], tol) && count[j] < itMax)
            fn = f(sol[j])
            der = fd(sol[j], h)
            der2 = f2d(sol[j], h)
            der3 = f3d(sol[j], h)
            if (der ≈ 0) && (der2 ≈ 0)
                # Third order Householder's method
                diff = (6*fn*der^2 - 3*fn^2*der2)/(6*der^3-6*fn*der*der2 + fn^2*der3)
            elseif (der ≈ 0)
                # Second order, Hailey's method
                diff = 2*fn*der/(-fn*der2 + 2*der^2)
            else
                # First order, Newton's method
                diff = fn/der
            end
            sol[j] -= diff
            count[j] += 1
        end
        if (count[j] == itMax)
            print("Maximum iterations exceeded and the amount by which ")
            println("Newton's last updated the solution was: ", diff)
        end
    end

    # Delete NaN entries from sol and corresponding count entries
    # While loop is used because sol's length will change over this loop
    j = 1;
    while j < length(sol)
        if isnan(sol[j])
            deleteat!(sol, j)
            deleteat!(count, j)
        end
        j += 1;
    end
    return sol, count
end

"""
    findRoot(f::Function, h::Float, tol::Float, itMax::Integer, a::Number, b::Number, N::Integer)

Uses the bisection method to get an initial guess of the root(s) of f on the 
domain [a, b] with N subdivisions, then uses Newton's method with a maximum of 
itMax iterations and a relative error tolerance of tol. h is the step size 
used to approximate the derivative. 
"""
function findRoot(f, h, tol, itMax, a, b, N)
    initGuess = bisection(f, N, a, b)
    sol, count = newtons(f, h, tol, itMax, initGuess);
    return sol, count
end

"""
    printSoln(sol::Vector{Float}, count::Vector{Integer})

Print the solution and the number of iterations used.
"""
function printSoln(sol, count)
    if (length(sol) == 1)
        println("Root = ", sol[1], ", count = ", count[1])
    else
        for i=1:length(sol)
            println("The $(i)th root = ", sol[i], ", count = ", count[i])
        end
    end
end

# This is where you specify the function you want to
# find the root of
function f(x)
    return x^4 + 3x^3 + 4x^2 + 3x + 1
end

# Initialize required variables
h = 1e-10
tol = 1e-15
itMax = 1000
a = -100
b = 100
N = 100000
sol, count = findRoot(f, h, tol, itMax, a, b, N)
printSoln(sol, count)