# Still Q14
x = 2;
i = 0;

function f(x)
    return tan(x) + 2*tanh(x)
end

function fd(x)
    return (sec(x))^2+2*(sech(x))^2
end

diffn = f(x)/fd(x);

while abs(diffn) > 1e-15
    diffn = f(x)/fd(x);
    x    += - diffn;
    i    += 1;
end

res=f(x);

println("$(i) iterations were used to solve this problem.")
println("x is $(x).")
println("f(x) is $(res).")