function [Zinp] = InputImpedanceOneVessel(NumModes,omegan,L,Z0,c,Gamma)
    % input impednace for Womersley solution
    % only for perturbation part
    % Zinp(1)=0 is for the steady part, zero frequency

    % in frequency domain
    Zinp = zeros(1,NumModes);
    for k=2:NumModes
        Zinp(k) = Z0(k) * (1.0 + Gamma(k)*exp(-1i*omegan(k)*2.0*L/c(k)))...
            /(1.0 - Gamma(k)*exp(-1i*omegan(k)*2.0*L/c(k)));  
    end
   
end