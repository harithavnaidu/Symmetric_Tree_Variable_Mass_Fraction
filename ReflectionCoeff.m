function [Gamma] = ReflectionCoeff(NumModes,ZT,Z0)
    %get reflection coefficient as defined
    
    Gamma = zeros(1,NumModes);
    for k = 2:NumModes
        if (ZT(k) + Z0(k)==0.0 )
            disp('exit: zero = ZT + Z0');
            break
        else
            Gamma(k) = (ZT(k) - Z0(k))/(ZT(k) + Z0(k));
        end
    end
end