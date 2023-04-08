function y=dWmdx_act(x)

global rho_w S lmax lmin

lm=lmax;
l0=lmin;

y=S/rho_w*(1-(lm-x)^2/(lm-l0)^2);

