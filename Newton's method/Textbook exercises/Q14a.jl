# Q14
using Plots;
x=LinRange(pi/2+0.1,3*pi/2-1e-1,10001);
y=tan.(x) + 2*tanh.(x);
plot(x,y)