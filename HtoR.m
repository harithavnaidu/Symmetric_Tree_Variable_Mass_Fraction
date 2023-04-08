function H = HtoR(R)

R = R*10^6;
H = 140/1000*R+1.26;

H = H/R;