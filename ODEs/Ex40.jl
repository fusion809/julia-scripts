# Ex 40
function f(t,x)
    return t*x
end

# Function that integrates the problem using RKF45
"""
	rkf45(h::Float64, epsilon::Float64, t0::Number, tf::Number, x0::Number)

uses the 4/5th order Runge-Kutta-Fehlberg method to integrate `f(t,x)` from t0 to tf with the initial condition
`x(0) = x0` and tolerance level of `epsilon`.
"""
function rkf45(h::Float64, epsilon::Float64, t0::Number, tf::Number, x0::Number)
    t = Float64[t0];
    x = Float64[x0];
    error_in_x = Float64[0];
    i = 1;
    while t[i]<tf
        h = min(h, tf-t[i]);
        k1 = h*f(t[i],x[i]);
        k2 = h*f(t[i]+h/4, x[i]+k1/4);
        k3 = h*f(t[i]+3*h/8, x[i]+3*k1/32+9*k2/32);
        k4 = h*f(t[i]+12*h/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197);
        k5 = h*f(t[i]+h, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104);
        k6 = h*f(t[i]+h/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40);
        x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        R = abs(x1-x2)/h;
        s = 0.84*(epsilon/R)^(1/4);
        if R<=epsilon
            push!(t, t[i]+h)
            push!(x, x1)
            i = i+1;
            push!(error_in_x, abs(x[i]-2*exp(1/2*t[i]^2)));
            h = s*h;
        else
            h = s*h;
        end
    end
    return [t, x, error_in_x]
end

h = 0.05;
t0 = 0;
tf = 0.25;
epsilon = 1e-14;
x0 = 2;
t = t0:h:tf;
N = length(t);
x_euler = zeros(N); x_euler[1] = x0;
x_rk4 = zeros(N); x_rk4[1] = x0;
x_rkf45 = 0;

for i = 1:N-1
    x_euler[i+1] = x_euler[i] + h*(x_euler[i]*t[i]);
    k1 = h*x_rk4[i]*t[i];
    k2 = h*(x_rk4[i]+1/2*k1)*(t[i]+1/2*h);
    k3 = h*(x_rk4[i]+1/2*k2)*(t[i]+1/2*h);
    k4 = h*(x_rk4[i]+k3)*(t[i]+h);
    x_rk4[i+1] = x_rk4[i] + 1/6*(k1+2*k2+2*k3+k4);
end

t_rkf45, x_rkf45, error_in_x_rkf45 = rkf45(h, epsilon, t0, tf, x0);

x_analytical = 2*exp.(1/2*t.^2);
error_in_euler = abs.(x_euler - x_analytical);
error_in_rk4 = abs.(x_rk4 - x_analytical);
rms_in_euler = sqrt(error_in_euler[2:end]'*error_in_euler[2:end]/(length(error_in_euler)-1));
rms_in_rk4 = sqrt(error_in_rk4[2:end]'*error_in_rk4[2:end]/(length(error_in_rk4)-1));
rms_in_rkf45 = sqrt(error_in_x_rkf45[2:end]'*error_in_x_rkf45[2:end]/(length(error_in_x_rkf45)-1));
["Data" "Euler's" "Runge-Kutta 4th order" "Runge-Kutta-Fehlberg"; 
"Root mean square of error" rms_in_euler rms_in_rk4 rms_in_rkf45; 
"Number of steps" length(x_euler)-1 length(x_rk4)-1 length(x_rkf45)-1]
