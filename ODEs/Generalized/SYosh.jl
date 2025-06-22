# StaticArrays help optimize this code
using StaticArrays;

function SYosh(fp::Function, fq::Function, params::NamedTuple, t0::Float64, tf::Float64, conds::SVector, dt::Float64)
    N = round(Int, (tf-t0)/dt, RoundUp);
    t = LinRange(t0, tf, N);
    h = step(t); 
    dims = round(Int, length(conds)/2);
    q = zeros(N, dims);
    p = zeros(N, dims);

    q[1,:] = conds[1:dims];
    p[1,:] = conds[(dims+1):length(conds)];
    # Coefficients
    w1 = 1 / (2 - 2^(1/3))
    w2 = -2^(1/3) / (2 - 2^(1/3))
    a = [w1/2, (w2 + w1)/2, (w2 + w1)/2, w1/2]
    b = [w1, w2, w1]

    for i in 1:N-1
        qi = q[i,:];
        pi = p[i,:];

        for j in 1:3
            pi += a[j] * h * fp(params, qi)
            qi += b[j] * h * fq(params, pi)
        end
        pi += a[4] * h * fp(params, qi)
        q[i+1,:] = qi;
        p[i+1,:] = pi; 
    end
    return t, q, p;
end