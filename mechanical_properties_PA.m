function KC = mechanical_properties_PA(p,R)

 global Ghe1 Ghe2 Ghc Ghm rho_w fiber_ang
 
% Passive wall
% elastin
Ghe1 = 1.27;% nondim.
Ghe2 = 1.27;% nondim.
Sigma_e = 28.8341;
Sigma_c = 178.5934;%p/(0.07); %120000;
Sigma_m = 24.5093;%15*133.32/(0.05/20*k+0.1);

c1 =  Sigma_e / ((Ghe1^2-(Ghe1^2*Ghe2^2)^(-1))*(0.3*rho_w));

c3 = 1.05; % for collagen fibers
c5 = 0.75;  % for passive smooth muscle

% fiber_ang = 45*pi/180; % angle of helical collagen fibers
Ghc = 1.154;% nondim. %HG changed from 1.034. June 22nd
str_ch = 0.0;
ang = [0.0 90.0 45.0 90+45.0]*(pi/180);
c_frac0 = [.2 .2 .3 .3];
for i=1:4
    str_ch = str_ch + c_frac0(i)*Ghc^2*(Ghc^2-1)*...
        exp(c3*(Ghc^2-1.0)^2)*(sin(ang(i)))^2;
end
c2 = Sigma_c/ (str_ch*(0.3*rho_w));

Ghm = 1.21;% nondim.
VSM_act = dWmdx_act(1)*0.3*rho_w;
c4 = (Sigma_m-VSM_act)/(Ghm^2*(Ghm^2-1.0)*exp(c5*(Ghm^2-1.0)^2)*(0.3*rho_w));

KC = [c1, c2, c3, c4, c5];

end