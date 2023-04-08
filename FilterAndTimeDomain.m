% Get function transformation from frequency to time domain: F -> f
% Use filtering of messy modes(need to check for each flow):
%   - trancation of NumModes frequency modes 
% Only for: N > 2*NumModes-2
% Note: may correct filtering to remove zero modes in between nozero modes!
%
% (by Vasilina 06-19-2018)
function [fTotal]=FilterAndTimeDomain(NumModes,N,Fn)

    F = zeros(N,1);
    F(1:NumModes) = Fn;
    
    % symmetric conjugate
    % F(k) = F(N-k+2), k=2,...N, F(2)=F(N), F(10)=F(N-8)
    F(N-NumModes+2:N) = flip(conj(Fn(2:NumModes)));

    % fTotal = N*ifft(F); 
    fTotal = N*ifft(F,'symmetric'); 

end