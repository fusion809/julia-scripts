# Chebyshev decomposition of a function on [t0, t1]
t0=1;
t1=10;
N=1000;
n=0:1:N;
x=-cos.(pi*n/N);
T=cos.(acos.(x)*n');
t=((t1-t0)/2)*x.+(t0+t1)/2;
f=exp.(cos.(t));
a=T\f;
