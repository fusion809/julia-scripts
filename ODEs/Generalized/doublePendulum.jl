using PyPlot;
include("RKF45.jl");

function DP(params::NamedTuple, t::Float64, vars::SVector{4, Float64})::SVector{4, Float64}
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
params = (g = 9.81, l1 = 1.0, l2 = 1.0, m1 = 1.0, m2 = 1.0);

# Initial conditions and domain of integration
t0 = 0.0;
tf = 60.0;
theta10, ptheta10, theta20, ptheta20 = pi/2, 0.0, pi/2, 0.0;
conds = @SVector [theta10, ptheta10, theta20, ptheta20];

# Error tolerance and initial step size
epsilon = 1e-10;
dtInitial = 0.1;

# Solve problem and extract solution values
@time begin
solution = RKF45(DP, params, t0, tf, conds, epsilon, dtInitial);
end
t, vars = solution.t, solution.vars;
theta1, ptheta1, theta2, ptheta2 = vars[:,1], vars[:,2], vars[:,3], vars[:,4];
# Positions of the pendulum bobs
x1 = l1*sin.(theta1);
y1 = -l1*cos.(theta1);
x2 = x1 + l2*sin.(theta2);
y2 = y1 - l2*cos.(theta2);

# Plots
PyPlot.figure(1);
PyPlot.clf();
PyPlot.plot(theta1, ptheta1, label=L"p_{\theta_1}");
PyPlot.xlabel(L"\theta_1");
PyPlot.ylabel(L"p_{\theta_1}");
PyPlot.legend();
PyPlot.figure(2);
PyPlot.clf();
PyPlot.plot(theta2, ptheta2, label=L"p_{\theta_2}");
PyPlot.xlabel(L"\theta_2");
PyPlot.ylabel(L"p_{\theta_2}");
PyPlot.legend();
PyPlot.figure(3);
PyPlot.clf();
PyPlot.plot(t, vars[:,1], label=L"\theta_1");
PyPlot.plot(t, vars[:,2], label=L"p_{\theta_1}");
PyPlot.plot(t, vars[:,3], label=L"\theta_2");
PyPlot.plot(t, vars[:,4], label=L"p_{\theta_2}");
PyPlot.xlabel(L"t")
PyPlot.legend()
# Bob #1
PyPlot.figure(4);
PyPlot.clf();
PyPlot.plot(x1, y1, label="Pendulum bob 1 location")
PyPlot.legend()
# Bob #2
PyPlot.figure(5);
PyPlot.clf();
PyPlot.plot(x2, y2, label="Pendulum bob 2 location")
PyPlot.legend()