function p = tauv(d)

n = length(d);

v = get_vn(n);

p = -0.5*v'*d*v;