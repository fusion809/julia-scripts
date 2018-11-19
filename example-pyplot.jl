using Pkg;
Pkg.add("PyPlot")
using PyPlot;
n=0:1:10 ;
x=cos.(pi*n/10);
y=acos.(x*1/10);
PyPlot.plot(x,y)
