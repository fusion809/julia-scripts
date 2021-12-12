#!/usr/bin/env julia
N = 10000
x = LinRange(-5, 0, N+1)
tol = 1e-14
itMax = 100
roots = [-sqrt(5)/2-3/2, sqrt(5)/2-3/2]
change = zeros(size(x))

function f(x)
    return x^2+3*x+1
end

function fd(x)
    return 2*x+3
end

fval = f.(x)

for i=2:N+1
    if (sign(fval[i]) + sign(fval[i-1]) == 0)
        change[i] = 1
    end
end

initGuess = x[change.==1]
sol = initGuess
println("Our initial guesses = ", sol)
for j=1:length(initGuess)
    count = 0
    while (minimum(abs.(-roots.+sol[j])) > tol && count < itMax)
        global sol[j] -= f(sol[j])/fd(sol[j])
        println("x = ", sol[j])
        count += 1
    end
end
println("Exact roots are = ", roots)