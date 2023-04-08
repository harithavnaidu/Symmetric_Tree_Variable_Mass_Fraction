function [ZT] = TerminalImpedanceVessel(NumModes,Z0,Gamma)
    %  - counterpart function of ReflectionCoeff -
    % get terminal impedance from reflection coefficient Gamma
    % and characteristic impedance Z0
    
    ZT = zeros(1,NumModes);
    for k=2:NumModes
        if (Gamma == 1.0)
            ZT(k)= NaN;
        dips('Warning: vessel has a closed end, ifinite terminal impedance')
        else
            ZT(k) = Z0(k)*(1.0 + Gamma(k))/(1.0 - Gamma(k));
        end
    end
    
%     %trancation of 10 frequency modes
%     ZT_k = zeros(1,Nt);
%     for k=1:NumModes-1
%         ZT_k (1,k+1) = ZT(k+1);
%         ZT_k (1,Nt-k+1) = conj(ZT(k+1));
%     end
% 
%     %time domain
%     zTt = ifft(ZT_k);
    
end