# Q15 continued
x = 3.7;
i = 0;

function f(x)
    return x^4 - 4*x^3 + x^2 + 1.2
end

function fd(x)
    return 4*x^3 - 12*x^2 + 2*x
end

diffn = f(x)/fd(x);

while abs(diffn) > 1e-15
    diffn = f(x)/fd(x);
    x    += -diffn;
    i    += 1;
end

println("$(i) iterations were used to solve this problem.")
println("x is $(x).")
println("f(x) is $(res).")
