# Q13
x=2.5;
i=0;

function f(x)
    return 23-x^3
end

function fd(x)
    return -3*x^2
end

diffn=f(x)/fd(x);

while abs(diffn) > 1e-15
    diffn = f(x)/fd(x);
    x    += - diffn;
end

res = f(x);

println("$(i) iterations were used to solve this problem.")
println("x is $(x).")
println("f(x) is $(res).")