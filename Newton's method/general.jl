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
    change = zeros(size(x))
    fval = f.(x)
    for i=2:N+1
        if (sign(fval[i]) + sign(fval[i-1]) == 0)
            change[i] = 1
        end
    end

    initGuess = x[change.==1]
    return initGuess
end

"""
    newtons(f::Function, h::Float, tol::Float, itMax::Integer, 
    initGuess::Vector{Float})

Uses Newton's method to refine initGuess, our initial guess of the root(s) of 
f, until either itMax iterations has been performed or the relative magnitude 
of the update Newton's provides to the estimate is below tol. h is the step 
size used to approximate the derivative.  
"""
function newtons(f, h, tol, itMax, initGuess)
    function fd(x, h)
        return (f(x+h)-f(x-h))/(2*h)
    end
    sol = initGuess
    count = Int64.(zeros(size(initGuess)))
    for j=1:length(initGuess)
        diff = 1
        while (abs(diff/sol[j]) > tol && count[j] < itMax)
            diff = f(sol[j])/fd(sol[j], h)
            sol[j] -= diff
            count[j] += 1
        end
        if (count[j] == itMax)
            print("Maximum iterations exceeded and the amount by which ")
            println("Newton's last updated the solution was: ", diff)
        end
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

# This is where you specify the function you want to
# find the root of
function f(x)
    return x^4 + x^3 - 10x^2 - 4x + 16
end

h = 1e-10
tol = 1e-15
itMax = 1000
a = -100
b = 100
N = 100000
sol, count = findRoot(f, h, tol, itMax, a, b, N)
for i=1:length(sol)
    println("The $(i)th root = ", sol[i], " count = ", count[i])
end