function r_int=hubradius(m,N,h_a_star,c_star)
r=N*m/2;
r_d=r-(h_a_star+c_star)*m;
if N>24
    N=24;
    r_int=r*(0.5+0.0344*sqrt(N-12));
else
    r_int=r*(0.5+0.0344*sqrt(N-12));
end
end