function f(t, r, rdot)
    G = 6.67408e-11;
    M = 5.97237e24;
    average_theta_dot = 2*pi/(27.321661);
    average_r = 384399e3;
    massless_momentum = (average_r^2)*average_theta_dot;
    return [rdot, (massless_momentum^2)/r^3 - G*M/(r^2)];
end

function rkf45(h, epsilon, t0, tf, r0, rdot0)
    r = Float64[r0];
    rdot = Float64[rdot0];
    t = Float64[t0];
    i = 1;

    while t[i] < tf
        h = min(h, tf-t[i]);
        K1 = h*f(t[i],r[i], rdot[i]);
        k1 = K1[1];
        l1 = K1[2];
        K2 = h*f(t[i]+h/4, r[i]+k1/4, rdot[i]+l1/4);
        k2 = K2[1];
        l2 = K2[2];
        K3 = h*f(t[i]+3*h/8, r[i]+3*k1/32+9*k2/32, rdot[i]+3*l1/32+9*l2/32);
        k3 = K3[1];
        l3 = K3[2];
        K4 = h*f(t[i]+12*h/13, r[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, rdot[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197);
        k4 = K4[1];
        l4 = K4[2];
        K5 = h*f(t[i]+h, r[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, rdot[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104);
        k5 = K5[1];
        l5 = K5[2];
        K6 = h*f(t[i]+h/2, r[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, rdot[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40);
        k6 = K6[1];
        l6 = K6[2];
        r1 = r[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        rdot1 = rdot[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        r2 = r[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        rdot2 = rdot[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
        R = abs(r1-r2)/h;
        s = 0.84*(epsilon/R)^(1/4);
        if R<=epsilon
            push!(t, t[i]+h)
            push!(r, r1)
            push!(rdot, rdot1)
            i = i+1;
            h = s*h;
        else
            h = s*h;
        end
    end
    return [t, r, rdot];
end

t0 = 0;
tf = 100;
r0 = 384399e3;
rdot0 = (405400e3-362600e3)/(27.321661);
epsilon = 1e-5;
h = 0.01;
t, r, rdot = rkf45(h, epsilon, t0, tf, r0, rdot0);
