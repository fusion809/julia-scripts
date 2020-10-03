using PyPlot;
include("RKF45.jl");

struct paramObj
    g::Float64
    l1::Float64
    l2::Float64
    m1::Float64
    m2::Float64
end

function f(params::paramObj, t, vars)
    g = params.g;
    l1 = params.l1;
    l2 = params.l2;
    m1 = params.m1;
    m2 = params.m2;
    theta1 = vars[1];
    ptheta1 = vars[2];
    theta2 = vars[3];
    ptheta2 = vars[4];

    # Variables to simplify later definitions
    C1 = ptheta1*ptheta2*sin(theta1-theta2)/(l1*l2*(m1+m2*(sin(theta1-theta2)^2)));
    C2 = ((l2^2)*m2*ptheta1^2 + (l1^2)*(m1+m2)*ptheta2^2 - l1*l2*m2*ptheta1*ptheta2*cos(theta1-theta2))/(2*(l1^2)*(l2^2)*(m1+m2*sin(theta1-theta2)^2))*sin(2*(theta1-theta2));
    
    # Derivatives
    thetaDot1 = (l2 * ptheta1 - l1*ptheta2 * cos(theta1))/((l1^2)*l2*(m1+m2*((sin(theta1-theta2))^2)));
    pthetaDot1 = -(m1+m2)*g*l1*sin(theta1) - C1 + C2;
    thetaDot2 = (l1*(m1+m2)*ptheta2-l2*m2*ptheta1*cos(theta1-theta2))/(l1*(l2^2)*m2*(m1+m2*(sin(theta1-theta2)^2)));
    pthetaDot2 = -m2*g*l2*sin(theta2)+C1-C2;

    # Return statement
    return [thetaDot1, pthetaDot1, thetaDot2, pthetaDot2];
end

# Problem parameters
g = 9.81;
l1 = 1.0;
l2 = 1.0;
m1 = 1.0;
m2 = 1.0;
params = paramObj(g, l1, l2, m1, m2);

# Initial conditions and domain of integration
t0 = 0.0;
tf = 10.0;
theta10 = pi/2;
ptheta10 = 0;
theta20 = pi/2;
ptheta20 = 0;
conds = [theta10 ptheta10 theta20 ptheta20];

# Error tolerance and initial step size
epsilon = 1e-10;
dtInitial = 0.1;

# Solve problem and extract solution values
solution = RKF45(f, params, t0, tf, conds, epsilon, dtInitial);
vars = solution.vars;
t = solution.t;
theta1 = vars[:,1];
ptheta1 = vars[:,2];
theta2 = vars[:,3];
ptheta2 = vars[:,4];

# Plots
PyPlot.figure(1);
PyPlot.clf();
PyPlot.plot(theta1, ptheta1);
PyPlot.figure(2);
PyPlot.clf();
PyPlot.plot(theta2, ptheta2);
PyPlot.figure(3);
PyPlot.clf();
PyPlot.plot(t, vars);