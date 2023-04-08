% updated 10-17-17
%
% get fequency-dependent and complex values of 
% wave speed and M,g coefficients        
% from deformable Womersley solution for oscillations
function [cn,Mn,gn,c_Rn] = WomersleySolutionCoeff(NumModes,omegan,mu,rho,... 
                                             R,E,sigma,h,rho_s)
   cn = zeros(NumModes,1);
   Mn = zeros(NumModes,1);
   gn = zeros(NumModes,1);
   
   % first component is dummy

    % Inviscid wave speed
    c_0 = sqrt(E*h/(2*rho*R));

    %---------frequency, Womersley #, Lambda and g factor vectors----------
    alphan = zeros(NumModes,1);
    Lambdan = zeros(NumModes,1);
    for k=2:NumModes;
%         omegan(k) = 2.0*pi*(k-1)/T;
        alphan(k) = sqrt(rho*omegan(k)/mu)*R;
        Lambdan(k,1) = alphan(k,1)*1i^1.5;
        gn(k) = 2.0*besselj(1,Lambdan(k))/(Lambdan(k)*besselj(0,Lambdan(k)));
    end
    
    %------------- Womersley deformable wall solution----------------------
    % coefficients of the complex polynomial equation
    % quadratic equation, 3 coefficients, 2 roots
    coef1n = zeros(NumModes,1);
    coef2n = zeros(NumModes,1);
    coef3n = zeros(NumModes,1);
    for k=2:NumModes;
        coef1n(k) = (gn(k)-1.0)*(sigma^2-1.0);
        coef2n(k) = rho_s*h*(gn(k)-1.0)/(rho*R)+(2.0*sigma-0.5)*gn(k)-2.0;
        coef3n(k) = 2.0*rho_s*h/(rho*R)+gn(k);
    end
    freq_coeffsn = zeros(NumModes,3);
    for k=2:NumModes;
        freq_coeffsn(k,:) = [coef1n(k) coef2n(k) coef3n(k)];
    end
    vin = zeros(NumModes,2);
    vin_1 = zeros(NumModes,1);
    vin_2 = zeros(NumModes,1);
    for k=2:NumModes
        vin(k,:) = roots(freq_coeffsn(k,:));
        vin_1(k) = vin(k,1);
        vin_2(k) = vin(k,2);
    end
    
    % ------------------complex wave speeds -------------------------------
    c_coeff_1n = zeros(NumModes,1);
    c_coeff_2n = zeros(NumModes,1);
    for k=2:NumModes
        c_coeff_1n(k) = sqrt(2.0/((1.0-sigma*sigma)*vin_1(k)));
        c_coeff_2n(k) = sqrt(2.0/((1.0-sigma*sigma)*vin_2(k)));
    end
    
    % wave speeds
    c_1n = zeros(NumModes,1);
    c_2n = zeros(NumModes,1);
    for k=2:NumModes
        c_1n(k) = c_coeff_1n(k)*c_0;
        c_2n(k) = c_coeff_2n(k)*c_0;
    end
    
    % |c_1| < c_0, |c_2| >c_0
    % c_1 is the right solution, since it is the one whose modulus is < 1. 
    % This is a physical condition, since the wave speed has to be lower
    % than that in an inviscid fluid
    vinc = zeros(NumModes,1);
    for k=2:NumModes
        if (abs(c_1n(k)) <= c_0)
            cn(k) = c_1n(k);
            vinc(k) = vin_1(k);
        elseif ((abs(c_2n(k)) <= c_0) && (abs(c_2n(k)) <= abs(c_1n(k))))  %??
    %     elseif (le(abs(c_2n),c_0) && le(abs(c_2n),abs(c_1n))) 
            cn(k) = c_2n(k);
            vinc(k) = vin_2(k);
        end
    end
    
    % --------------- Estimate of the spatial wavelength-------------------
    inverse_cn = zeros(NumModes-1,1);
    c_Rn = zeros(NumModes,1);
    c_In = zeros(NumModes,1);
    for k=2:NumModes;
        inverse_cn(k) = cn(k)^(-1);
        c_Rn(k) = 1.0/real(inverse_cn(k));
        c_In(k) = 1.0/imag(inverse_cn(k));
    end
%     Lwave = c_Rn(1)*T;
    
    % -------------------------get M --------------------------------------
    for k=2:NumModes
       Mn(k) = (2.0 + vinc(k)*(2.0*sigma-1.0))/(vinc(k)*(2.0*sigma-gn(k)));
    end
end