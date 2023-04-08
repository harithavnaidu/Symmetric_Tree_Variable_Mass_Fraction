function y=dWmdx(x,kc)

% global kc

y=kc(4)*(x*(x^2-1))*exp(kc(5)*(x^2-1)^2);

