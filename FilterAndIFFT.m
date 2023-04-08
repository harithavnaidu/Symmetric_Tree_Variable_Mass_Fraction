% Get function transformation from frequency to time domain: F -> f
% Use filtering of messy modes(need to check for each flow):
%   clean from zero frequencies of Qn, excluding steady, 
% Total = steady state + oscillatory contribution
% Oscillatory - without zero frequency mode (k=1)
%   %trancation of 10 frequency modes 
% only for: N > 2*NumModes-2
% Note: need to correct filtering to remove zero modes in between nozero modes!!!

function [fTotal,fOscil]=FilterAndIFFT(NumModes,N,Fn)

    F = zeros(N,1);
    F(1:NumModes) = Fn;
    
    % symmetric conjugate
    % F(k) = F(N-k+2), k=2,...N, F(2)=F(N), F(10)=F(N-8)
    F(N-NumModes+2:N)=flip(conj(Fn(2:NumModes)));

    %     fTotal = N*ifft(F); 
    fTotal = N*ifft(F,'symmetric'); 
    fOscil = fTotal - F(1);

end