% modified by Vasilina on August 2, 2018
function [q, p_term] = SymmetricArterialTree2(N_gen,q_parent,p_terminal,Res)
% Find steady flow and terminal pressure at each generation
% symmetric bifurcation: half split of flow
% p(k) = p(k+1) + dp(k+1); dp(k)=q/(2^(k-1))*Res(k)

p_term = zeros(N_gen,1);
q = zeros(N_gen,1);

p_term(N_gen) = p_terminal;
q(N_gen) = q_parent/(2^(N_gen-1));

for k = N_gen-1:-1:1
    q(k) = q_parent/(2^(k-1));
    p_term(k) = p_term(k+1) + Res(k)*q(k);
end
