#!/usr/bin/env julia
using BigRationals

function newtons(f, params, x0, itMax, tol)
    x = [x0];
    count = 0;
    diff = 1;
    dx = 1-10;

    function fd(f, params, x, dx)
        dfdx = (f(x+dx, params)-f(x-dx, params))/(2*dx)
        return dfdx
    end

    while (abs(diff) > tol && count < itMax)
        diff = f(x[end], params)/fd(f, params, x[end], dx)
        x = push!(x, x[end] - diff)
        count += 1
    end
    return x
end

function sqrtZero(x, a)
    func = (x^2-a)
    return func
end

function sqrtApp(a, itMax, tol)
    a = BigRational(a)
    upPerfSqrt = BigRational(ceil(sqrt(a)))
    btPerfSqrt = BigRational(floor(sqrt(a)))
    upPerfSq = upPerfSqrt^2
    btPerfSq = btPerfSqrt^2
    if (upPerfSq - a) >= (a - btPerfSq)
        xTang = btPerfSqrt + (a - btPerfSq)/(2*btPerfSqrt)
    else
        xTang = upPerfSqrt + (a - upPerfSq)/(2*upPerfSqrt)
    end
    x = newtons(sqrtZero, a, xTang, itMax, tol)
    return x
end
