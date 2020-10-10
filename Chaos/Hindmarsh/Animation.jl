# This is largely copied from code in the docs:
# http://docs.juliaplots.org/latest/ (Author: Thomas Breloff)
# My attempt at animation failed see Lorenz.jl and 
# https://discourse.julialang.org/t/how-to-create-an-animation-based-on-t-x-y-and-z-data/25672/4
# define the Hindmarsh-Rose attractor
using Plots;
mutable struct Hindmarsh
    dt; a; b; c; d; r; s; xr; I; x; y; z
end

function step!(l::Hindmarsh)
    dx = l.y-l.a*(l.x)^3+l.b*(l.x)^2-l.z+l.I       ; l.x += l.dt * dx
    dy = l.c-l.d*(l.x)^3-l.y                       ; l.y += l.dt * dy
    dz = l.r*(l.s*(l.x-l.xr)-l.z)                  ; l.z += l.dt * dz
end

attractor = Hindmarsh((dt = 0.001, a=1., b=3., c=1., d=5., r=1e-3, s=4., xr=-1.6, I=-2.0, x = 1., y = 1., z = 1.)...)


# initialize a 3D plot with 1 empty series
plt = plot3d(1, xlim=(-25,25), ylim=(-25,25), zlim=(0,50),
                title = "Hindmarsh-Rose Attractor", marker = 2, 
                xlabel="x", ylabel="y", zlabel="z",
                size=(1920,1080),
                label="")

# build an animated gif by pushing new points to the plot, saving every 10th frame
@gif for i=1:10000
    step!(attractor)
    push!(plt, attractor.x, attractor.y, attractor.z)
end every 100