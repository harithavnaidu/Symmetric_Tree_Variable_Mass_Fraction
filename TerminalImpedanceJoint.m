function [ZT] =  TerminalImpedanceJoint(NumModes,Zinp1,Zinp2)
    % get terminal impedance of parent vessel
    % computed at joint
    
    ZT = zeros(1,NumModes);
    for k = 2:NumModes
        if ((Zinp1(k) + Zinp2(k))==0.0 )
            disp('exit: zero = Zinp1+Zinp2');
            break
        else
            ZT(k)= Zinp1(k)*Zinp2(k)/(Zinp1(k) + Zinp2(k));  
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