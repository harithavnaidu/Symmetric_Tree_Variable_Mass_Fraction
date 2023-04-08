function [Pn] = Precovery(NumModes,omegan,L,x,Hfn,cn,Gamma)
    
    %frequency domain   
    Pn = zeros(NumModes,1);
    for k=2:NumModes
        Pn(k) = Hfn(k) *exp(-1i*omegan(k)*x/cn(k)) * ...
            (1.0 + Gamma(k)*exp(1i*omegan(k)*2.0*(x-L)/cn(k)));
    end
   
end