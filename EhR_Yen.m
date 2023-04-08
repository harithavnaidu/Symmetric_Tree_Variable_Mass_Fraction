%% This script plots the data from Yen RT, Rong Z, Zhang B (1990) Elasticity of pulmonary blood vessels in human lungs.
EhrY_mean = 62.5*133.32; %Yen1990 (from Qureshi 2014) lambda=0.012 mmHg
for kk=1:Newgen
    if ( Radius(kk) <= 300/10^6 )
        lambda = 1.112/100/0.73555912101486/133.32; % unit: %/cmH2O*cmH2OtommHg
        EhrY(kk) = 3/(4*lambda);
    elseif (Radius(kk) > 300/10^6 && Radius(kk) <= 400/10^6)
        lambda = 0.913/100/0.73555912101486/133.32; % unit: %/cmH2O*cmH2OtommHg
        EhrY(kk) = 3/(4*lambda);
    elseif (Radius(kk) > 400/10^6 && Radius(kk) <= 600/10^6)
        lambda = 0.859/100/0.73555912101486/133.32; % unit: %/cmH2O*cmH2OtommHg
        EhrY(kk) = 3/(4*lambda);
    elseif (Radius(kk) > 600/10^6 && Radius(kk) <= 1000/10^6)
        lambda = 1.213/100/0.73555912101486/133.32; % unit: %/cmH2O*cmH2OtommHgz
        EhrY(kk) = 3/(4*lambda);
    else
        lambda = 1.255/100/0.73555912101486/133.32; % unit: %/cmH2O*cmH2OtommHg*mmHgtokPa
        EhrY(kk) = 3/(4*lambda);
    end
end