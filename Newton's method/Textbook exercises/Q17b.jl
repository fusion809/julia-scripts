# Q17 continued
# Root appears to be at 0.3
x=0.3;
i=0;

function f(x)
    return 3*x-exp(-x)
end

function fd(x)
    return 3+exp(x)
end

diffn = f(x)/fd(x);

while abs(diffn) > 1e-15
    x    += -diffn;
    diffn = f(x)/fd(x);
    i    += 1
end

println("$(i) iterations were used to solve this problem.")
println("x is $(x).")
println("f(x) is $(res).")