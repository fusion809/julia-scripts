r1_0        = 147098450e3
r1dot_0     = 0;
r2_0        = 362600e3
r2dot_0     = 0;
theta1_0    = 0;
thetadot1_0 = -(29782.7)/(r1_0);
theta2_0    = 0;
thetadot2_0 = -(1022)/(r2_0);
x2_0        = r1_0 * cos(theta1_0);
xdot2_0     = r1dot_0 * cos(theta1_0) - r1_0 * thetadot1_0 * sin(theta1_0);
y2_0        = r1_0 * sin(theta1_0);
ydot2_0     = r1dot_0 * sin(theta1_0) + r1_0 * thetadot1_0 * cos(theta1_0);
x3_0        = r2_0 * cos(theta2_0);
xdot3_0     = r2dot_0 * cos(theta2_0) - r2_0 * thetadot2_0 * sin(theta2_0);
y3_0        = r2_0 * sin(theta2_0);
ydot3_0     = r2dot_0 * sin(theta2_0) + r2_0 * thetadot2_0 * cos(theta2_0);
m1          = 5.972168e24;
m2          = 1.9885e30;
m3          = 7.346e22;
px2_0       = m2*xdot2_0;
py2_0       = m2*ydot2_0;
px3_0       = m3*xdot3_0;
py3_0       = m3*ydot3_0;

function fq(params, p)
    px2 = p[1];
    py2 = p[2];
    px3 = p[3];
    py3 = p[4];
    m2 = params.m2;
    m3 = params.m3;
    xdot2 = px2/m2;
    ydot2 = py2/m2;
    xdot3 = px3/m3;
    ydot3 = py3/m3;
    return [xdot2, ydot2, xdot3, ydot3];
end

function fp(params, q)
    x2 = q[1];
    y2 = q[2];
    x3 = q[3];
    y3 = q[4];
    G = 6.674e-11;
    m2 = params.m2;
    m3 = params.m3;
    r22 = (x2^2+y2^2)^(3/2);
    r32 = ((x3-x2)^2+(y3-y2)^2)^(3/2);
    px2dot = -G*m2*(m1*x2/r22 + m3*(x2-x3)/r32);
    py2dot = -G*m2*(m1*y2/r22 + m3*(y2-y3)/r32);
    px3dot = -G*m3*(m1*x3/r22 + m3*(y3-y2)/r32);
    py3dot = -G*m3*(m1*y3/r22 + m3*(y3-y2)/r32);
    return [px2dot, py2dot, px3dot, py3dot];
end

params = (m1=m1, m2=m2, m3=m3);
conds = @SVector [x2_0, y2_0, x3_0, y3_0, px2_0, py2_0, px3_0, py3_0];
dt = 1.0;
t0 = 0.0;
tf = 60*60*24*365.25;
include("SYosh.jl")
t1, x, p = SYosh(fp, fq, params, t0, tf, conds, dt);
x2 = x[:,1];
y2 = x[:,2];
x3 = x[:,3];
y3 = x[:,4];
px2 = p[:,1];
py2 = p[:,2];
px3 = p[:,3];
py3 = p[:,4];
using PyPlot;
PyPlot.figure(1)
PyPlot.plot(x2, y2, color=:red)
PyPlot.plot(x3, y3, color=:blue)