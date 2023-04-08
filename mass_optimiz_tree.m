% modified by Vasilina on August 2, 2018
function [Mtotal, lhat, phi_e, SMCtoCOL] = mass_optimiz_tree(p,Q,k)

global rho_w beta gamma R0 alpha_t Ghe1 Ghe2 Ghc Ghm  fiber_ang   

%% Initialization
lhat = 0.9^(k-1);  % why?

LamZ = 1;

%% Newton-Raphson optimization
err = 1000;
Y = lhat;

c=1;
while err>=1e-8
    %%
    
    kc = mechanical_properties_PA_2(p,lhat*R0);
    
    lhat = Y;
    [phi_e, phi_c, phi_m] = mass_fracs(2*lhat*R0);
    
    SMCtoCOL = phi_m/phi_c;
    
    beta_t = 1/LamZ*(1/(1+SMCtoCOL)*(0.2+ 0.3*sin(fiber_ang)^2 + ...
            0.3*sin(-fiber_ang)^2)*Ghc*dWcdx(Ghc,kc) +...
            SMCtoCOL/(1+SMCtoCOL)*(Ghm*dWmdx(Ghm,kc)+dWmdx_act(1)));
    
    phi_t = 1 - phi_e;
    
    alpha_act = 0.00872;
    
    VSM_act = dWmdx_act(1)*0.3*rho_w;
    
    alpha_total = alpha_t + alpha_act*VSM_act;
    alpha_act = 0;
        
    beta_e = 1/LamZ*dWedx(1, Ghe1, Ghe2,kc);
    %    active-tone cost +2*pi*alpha_act*SMCtoCOL/(SMCtoCOL+1)*phi_t*dWmdx_act(Ghm)/(phi_e*beta_e+phi_t*beta_t)*(3*p*lhat^2*R0^3)+...
    DF = 2*pi/(0.3*rho_w)*alpha_total*SMCtoCOL/(SMCtoCOL+1)*phi_t/(phi_e*beta_e+phi_t*beta_t)*(2*p*lhat*R0^2)+ ...
        2*pi/(0.3*rho_w)*alpha_t*1/(SMCtoCOL+1)*phi_t/(phi_e*beta_e+phi_t*beta_t)*(2*p*lhat*R0^2)+...
        alpha_act*SMCtoCOL/(SMCtoCOL+1)*phi_t*dWmdx_act(Ghm)/(0.3*rho_w) + ...
        2*beta*R0^2*lhat - 4*gamma*Q^2*R0^(-4)*lhat^(-5);
    
    D2F = 2*pi/(0.3*rho_w)*alpha_total*SMCtoCOL/(SMCtoCOL+1)*phi_t/(phi_e*beta_e+phi_t*beta_t)*(2*p*R0^2)+...
        2*pi/(0.3*rho_w)*alpha_t*1/(SMCtoCOL+1)*phi_t/(phi_e*beta_e+phi_t*beta_t)*(2*p*R0^2)+...
        2*pi*alpha_act*SMCtoCOL/(SMCtoCOL+1)*phi_t*dWmdx_act(Ghm)/(phi_e*beta_e+phi_t*beta_t)*(6*p*lhat*R0^3)+...
        2*beta*R0^2 + 20*gamma*Q^2*R0^(-4)*lhat^(-6);
    
    sol = DF/D2F;
    err = sqrt( sum(sol.^2)/sum(Y.^2) );
    Y = Y - sol;
    c=c+1;
end

lhat = Y;


[phi_e, phi_c, phi_m] = mass_fracs(2*lhat*R0);
SMCtoCOL = phi_m/phi_c;
beta_t = 1/LamZ*(1/(1+SMCtoCOL)*(0.2+ 0.3*sin(fiber_ang)^2 + 0.3*sin(-fiber_ang)^2)*Ghc*dWcdx(Ghc,kc) + SMCtoCOL/(1+SMCtoCOL)*(Ghm*dWmdx(Ghm,kc)+dWmdx_act(1)));
phi_t = 1 - phi_e;

Mtotal = (p*R0*lhat)/(phi_e*beta_e+phi_t*beta_t);

disp(['     Mass-Radius Optimization at gen. ',num2str(k),': Newton-Raphson iterations # ',num2str(c),', tollerance 1e-8']);
