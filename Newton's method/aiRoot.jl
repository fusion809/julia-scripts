using SpecialFunctions;
x = -2.338;
i = 0;
diffn = airyai(x)/airyaiprime(x);

while abs(diffn) > 1e-15
    diffn = airyai(x)/airyaiprime(x);
    x    += -diffn;
    i    += 1;
end

println("$(i) iterations were used.")
println("x is $(x).")
println("diffn is $(diffn).")