ha = 0.1;
hb = 0.05;
hc = ha/10;
x0 = 2;
t0 = 1;
tf = 2;
ta = t0:ha:tf;
tb = t0:hb:tf;
tc = t0:hc:tf;
xa = zeros(length(ta));
xa[1] = 2;
xb = zeros(length(tb));
xb[1] = 2;
xc = zeros(length(tc));
xc[1] = 2;

for i=1:length(ta)-1
    xa[i+1] = xa[i] + ha*(xa[i]*ta[i]/(ta[i]^2+2));
end

for i=1:length(tb)-1
    xb[i+1] = xb[i] + hb*(xb[i]*tb[i]/(tb[i]^2+2));
end

for i=1:length(tc)-1
    xc[i+1] = xc[i] + hc*(xc[i]*tc[i]/(tc[i]^2+2));
end
