function [Zn] = CharacteristicImpedanceOneVessel(NumModes,R,rho,cn,Mn,gn)
    %compute characteristic impedance
    
    %frequency domain
    Zn = zeros(1,NumModes); 
%     for k = 1:NumModes-1
%         Zn(k+1)=rho*cn(k)/(pi*R^2*(1.0-Mn(k)*gn(k))); %complex values
%     end
    for k = 2:NumModes
        Zn(k)=rho*cn(k)/(pi*R^2*(1.0-Mn(k)*gn(k))); %complex values
    end

%     %trancation of frequency modes 
%     Zk = zeros(1,Nt);
%     for k=1:NumModes-1
%         Zk (1,k+1) = Zn(k+1);
%         Zk (1,Nt-k+1) = conj(Zn(k+1));
%     end
%     zt = ifft(Zk);  %is real value

end