# Q16 continued
# Smallest root appears to be ~-0.3
x=-0.3;
i=0;

function f(x)
    return x^4-3*x^3+2*x^2-3*x-1.6
end

function fd(x)
    return 4*x^3-9*x^2+4*x-3
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