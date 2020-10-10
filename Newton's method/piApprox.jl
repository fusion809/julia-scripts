# Use Newton's method to approximate pi by finding the zero of
# sin(x) in the vicinity of x=22/7
x = 22.0/7.0;

difference = sin(x)/cos(x);
i = 0;
while abs(difference) > 1e-15
    difference = sin(x)/cos(x);
    x += - difference;
    i += 1;
end

println("x is $(x)")
println("$(i) iterations were used")
println("difference is $(difference)")