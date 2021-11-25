#!/usr/bin/env julia
tol = BigFloat(1e-50)
eps = BigFloat(1e-1)
count = 0
itMax = 1e3
x0 = BigFloat(2.8)
a = BigFloat(8)
x = x0

# Apply Newton's method
while (count < itMax && abs(eps) > tol)
	global eps = (x^2-a)/(2*x)
	global x -= eps
	global count += 1
end

# Print variables
println("x     = ", x)
println("count = ", count)
println("eps   = ", eps)
println("x^2-a = ", x^2-a)