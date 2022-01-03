#!/usr/bin/env julia
function findLog(N, tol, x, y)
    count = 1
    eps = 1
    while (count < N && abs(eps) > tol)
        eps = 2*(x-exp(y))/(exp(y)+x)
        y += eps
        count += 1
    end
    return y
end

function testLog(N, tol, x, y)
    logVal = findLog(N, tol, x, y)
    diff = logVal - log(x)
    println("diff = $diff")
end